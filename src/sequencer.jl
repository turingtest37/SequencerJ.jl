"""

Constant to use when specifying the Square Euclidean or L2 distance metric for a sequencing run.

    The L2 metric measures the sum of squares of the Euclidean or Manhattan-style 
    distance between two vectors of values. Because it operates only on values and not 
    on their underlying grid or axis, the L2 metric is less sensitive to the positions of  
    values along their grid than the EMD and Energy metrics.
    ``
        L2(x,y) = Σ ||x - y||²
    ``
See [`Distances.SqEuclidean`](@Ref)

"""
const L2 = Distances.SqEuclidean()
const KLD = Distances.KLDivergence()
const WASS1D = EMD()
const ENERGY = Energy()
const ALL_METRICS = (L2, WASS1D, KLD, ENERGY)

# Smallest allowed weight/distance (instead of 0)
const ϵ = 1e-6

# dictionaries to hold intermediate results, keyed by (metric, scale)
EOSeg = Dict{Tuple,Any}() # list of elongation and orderings (BFS,DFS), one per Segment

EOAlgScale = Dict{Tuple,Any}() # elongation and orderings for the cumulated weighted distances

"""
    sequence(A::VecOrMat{T}; scales=(1,4), metrics=ALL_METRICS, grid=nothing) where {T <: Real}

    Analyze the provided `m x n` matrix (or m vectors of vectors n) by applying one or more 1-dimensional statistical metrics to 
    each column of data, possibly after dividing the data into smaller row sets. Columns are compared pairwise for each combination 
    of metric and scale to create a n x n distance matrix

    The paper that describes the Sequencer algorithm and its applications can be found 
    on Arxiv: [https://arxiv.org/abs/2006.13948].
    ```
@misc{baron2020extracting,
    title={Extracting the main trend in a dataset: the Sequencer algorithm},
    author={Dalya Baron and Brice Ménard},
    year={2020},
    eprint={2006.13948},
    archivePrefix={arXiv},
    primaryClass={cs.LG},
    year=2020
}
```
"""
function sequence(A::VecOrMat{T}; scales=(1,4), metrics=ALL_METRICS, grid=nothing) where {T <: Real}

    @info "Sequencing data with shape $(size(A)) using metric(s) $(metrics) at scale(s) $(scales)..."

# We need to ensure that A is, or becomes, oriented so that
# observation data run in columns. Each column is a separate data series.
# Each row typically represents a distinct point in time or space where the
# observed data were captured.
#
    @assert all(sign.(A) .>= 0) && all(!any(isnan.(A)) && !any(isinf.(A))) "Input data cannot contain negative, NaN or infinite values."

    # take a vector of vectors and cat it into a matrix
    A = A isa Vector ? hcat(A...) : A

    # replace zeros and Infs with values that work for the math
    clamp!(A, eps(), typemax(eltype(A)))

    grid = _ensuregrid!(grid, A)

    MST_all = []
    η_all = Float64[]

    for k in metrics
        alg = k
        for l in scales

            # summary distance matrix for segments
            Dklms = zeros(size(A,2), size(A,2))
            # Elongation per chunk
            ηs = []
            # BFS orderings per chunk
            orderings = []

            # split the data and grid into chunks
            S, G = _splitnorm(A, grid, l)

            # Each m row in S contains n segments of data,
            # one for each of n data series
            @inbounds for i in eachindex(S, G)

                # m is a chunk of the input data matrix
                m = S[i]
                # local grid with the same dim 1 dimension as m
                localgrid = G[i] # a vector

                # Some algorithms can be performed on arbitrary real-valued grids
                if alg in (EMD, Energy)
                    alg = alg((localgrid, localgrid))
                end

                # All the heavy lifting happens here, in the distance calculations
                Dklm = abs.(pairwise(alg, m; dims = 2)) .+ ϵ

                # Find the minimum spanning tree of the distance matrix
                MSTklm = _mst(Dklm) #MST

                # Use closeness centrality to find a starting point for graph walks
                startidx = elong_start_index(MSTklm)

                # half len / half width + 1
                η = elongation(MSTklm, startidx)

                # weight the distance matrix by its elongation factor
                Dklm_e = η .* Dklm

                # ensure we stay on the playing field...
                clamp!(Dklm_e, eps(), typemax(eltype(Dklm_e)))
                
                # update the summary matrix
                Dklms .+= Dklm_e

                # do a breadth-first search on the segment/chunk MST to get a local optimal sequence
                bfso = bfs_tree(MSTklm, startidx)

                push!(ηs, η)
                push!(orderings, bfso)
            end

            if length(Dklms) < 1
                @warn("Unable to create distance matrices from the given data.")
                return nothing
            end

            # Do the same MST-based operations over the weighted per-segment results
            Dkl = Dklms / sum(ηs) # Weighted average over all chunks for <alg,s>
            MSTkl = _mst(Dkl)
            stidx = elong_start_index(MSTkl) # Least central point of averaged MST
            ηkl = elongation(MSTkl, stidx)
            BFSkl = bfs_tree(MSTkl, stidx)

            # Store these as intermediate results
            push!(MST_all, MSTkl) # All MSTs (1 per alg,s combo)
            push!(η_all, ηkl) # All elongations (1 per alg,s combo)
            global EOSeg[(alg,l)] = (ηs, orderings)
            global EOAlgScale[(alg,l)] = (ηkl, BFSkl)

            @info "$(k) at scale $(l): η = $(ηkl)"

        end
    end

    # Now use the accumulated MSTs and elongations to create a final N x N proximity matrix

    # nb of columns in A
    N = size(A,2)
    # Pw is a sparse proximity matrix to be filled with MST elongation-weighted edge
    # distances (which are actually weights as well)
    Pw = _weighted_p_matrix(N, MST_all, η_all)
    # Invert the proximity matrix to get a final distance matrix
    # w/ zeros on the diagonal
    D = inv.(Pw)
    # final minimum spanning tree for analysis
    MSTD = _mst(D) #MST
    stidx = elong_start_index(MSTD) # Least central point of averaged MST

    # final elongation
    ηD = elongation(MSTD, stidx)

    # final optimal data sequence
    BFSD = bfs_tree(MSTD, stidx)

    @info "Final average elongation: $(ηD)"

    return MSTD, ηD, BFSD
end

"Ensure the grid is compatible with the data. Create a grid if one was not provided."
function _ensuregrid!(grid, A)::AbstractVector
    if isnothing(grid)
        # create a sensible grid if one was not provided
        grid = float(collect(axes(A,1)))
    elseif !(eltype(grid) isa AbstractFloat)
        grid = float(grid)
    end
    @assert length(grid) == size(A,1) "Grid length must match dims 1 (row count) of A. Got grid:$(length(grid)) and A:$(size(A,1))"
    grid
end


"Create a N x N proximity matrix from weighted edges of the given graph(s).
Edge weights are provided as η coefficient(s) >= 1."
function _weighted_p_matrix(N, graphs, ηs)
    P = spzeros(N,N)

# fill the proximity matrix with elongation-weighted reciprocal edge weights, treated as distances.
    @inbounds for idx in eachindex(graphs, ηs)
        g = graphs[idx]
        W = LightGraphs.weights(g)
        η = ηs[idx]
        for e in edges(g)
            i,j = src(e),dst(e)
            d = W[i,j] # weight as distance
            # proximity is the inverse of distance
            P[i,j] = P[j,i] += η * 1.0 / d
        end
    end
    # All P[i] were elongated using an η. Now divide by Ση to get the elongation-weighted average.
    Pw = P / sum(ηs)
    # Put Inf on the 1,1 diagonal to form a true proximity matrix
    Pw .+= sparse(Matrix(Inf*I, size(Pw)...))

    return Pw
end

"Return the index of the least central node in the given graph."
elong_start_index(g) = last(findmin(closeness_centrality(g)))

"Return the half-length or mean of the given path lengths (distances)."
_halflen(paths) = mean(paths)

"Return the half-width or mean count of unique values in the given path lengths."
_halfwidth(paths) = mean(count(v->v==k, paths) for k in unique(paths)) / 2

"""`julia`

`elongation(g, startidx)`

Returns the ratio of the graph half-length (mean of path distances) over the 
half-width, defined as the mean count of shortest paths from the center node of a
minimum spanning tree over the graph.
"""
function elongation(g, startidx)
    spaths = dijkstra_shortest_paths(g, startidx)
    # remove those annoying zeros and Infs that ruin it for everyone else
    hlen = clamp(_halflen(spaths.dists), eps(), floatmax(Float32))
    hwid = clamp(_halfwidth(spaths.dists), eps(), floatmax(Float32))
    return hlen / hwid + 1
end

"Return a weighted graph of the minimum spanning tree of the given distance matrix,
calculated using the Kruskal algorithm."
function _mst(dm::AbstractMatrix)

    # Zeros and Infs mess up the calculations...
    clamp!(dm, ϵ, typemax(eltype(dm)))

    # By forcing an unsymmetric distance matrix (from KLD and others) into a Symmetric one
# we are losing 1/2 (the lower half) of the distances
# TODO Should we construct a separate graph for the lower half as well?
    g = issymmetric(dm) ? SimpleWeightedGraph(dm) : SimpleWeightedGraph(Symmetric(dm))

    # Minimum Spanning Tree using the Kruskal algorithm
    mst = kruskal_mst(g; minimize=false) # calculate the maximum tree instead of the minimized version
    @debug "kruskal MST" mst
    Ii = src.(mst)
    Ji = dst.(mst)
    Vi = weight.(mst)
    @debug "I J V" Ii Ji Vi

# Populate a distance matrix with the MST nodes and edge weights
    S = sparse(Ii, Ji, Vi, size(dm)...)
    @debug "big beautiful S symmetric?" S issymmetric(S)

    g = issymmetric(S) ? SimpleWeightedGraph(S) : SimpleWeightedGraph(Symmetric(S))
    return SimpleWeightedGraph(g)
end


"Split the given `m x n` matrix and grid into `scale` parts along the column dimension such
that the resulting matrices (chunks) are approximately of dimension `k x n` where `k ≊ m / scale`"
function _splitnorm(A::AbstractMatrix, grid::AbstractVector, scale::Int)
    @debug "scale" scale
    @debug "grid" grid
    @assert scale <= length(A[:,1]) "Scale ($(scale)) cannot be larger than the number of data elements ($(length(A[:,1])))."

    collen = size(A,1)
    @debug "collen" collen
    # break up each column of A into N chunks, where N = scale
    chunklen = cld(collen, scale)
    @debug "chunklen" chunklen

    N = cld(collen, chunklen) # Number of chunks
    @debug "N (number of bins)" N

    slices=[]
    grids=[]
    # range of chunk starting indices
    chunkrange = 1:chunklen:collen
    for i in chunkrange
        ii = i + chunklen - 1 # thank you, 1-based indexing
        # use the "end" keyword on the last slice to mop up the remainder
        ii = clamp(ii, 1, collen)
        S = A[i:ii, :]
        @debug "slice S" S
        G = grid[i:ii]
        # now normalize each column of the slice
        Σslice = sum(S, dims=1)
        push!(slices, S ./ Σslice)
        push!(grids, G)
    end
    return slices, grids
end


