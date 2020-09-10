"""
    SequencerResult

Type whose fields contain the results of a Sequencer run.
```@example
using SequencerJ #hide
A = rand(10,10);
r = sequence(A;);
order(r)
```

"""
struct SequencerResult
    EOSeg::Dict{Tuple,Any} # list of elongation and orderings (BFS,DFS), one per Segment    
    EOAlgScale::Dict{Tuple,Any} # elongation and orderings for the cumulated weighted distances
    D::AbstractMatrix # final distance matrix
    mst::LightGraphs.AbstractGraph # final mst
    η::Real # final elongation
    order::AbstractVector #final ordering from bfs
end

"Alias for elong"
η(r::SequencerResult) = elong(r)

"""

    D(r::SequencerResult)

Return the final distance matrix from the Sequencer run.
"""
D(r::SequencerResult) = r.D

"""

    mst(r::SequencerResult)

Return the final minimum spanning tree graph from a Sequencer run. 
"""
mst(r::SequencerResult) = r.mst

"""

    elong(r::SequencerResult)
    
Return the elongation coefficient of the final, weighted graph.
"""
elong(r::SequencerResult) = r.η

"""

    order(r::SequencerResult)

Return the result column indices, as determined by the Sequencer algorithm.
"""
order(r::SequencerResult) = r.order

"Sensibly display a SequencerResult object."
show(io::IO, s::SequencerResult) = write(io, "Sequencer Result: η = $(@sprintf("%.4g", elong(s))), order = $(order(s)) ")

"Convert index to weight"
i2weight(x::AbstractVector) =  1 .- x ./ (maximum(x) + eps())

"""
    sequence(A::VecOrMat{T}; 
        scales=nothing,
        metrics=ALL_METRICS,
        grid=nothing,
        weightrows=false,
        rowfn=i2weight,
        ) where {T <: Real}

Analyze the provided `m x n` matrix (or m vectors of vectors n) by applying one or more 1-dimensional statistical metrics to 
each column of data, possibly after dividing the data into smaller row sets. Columns are compared pairwise for each combination 
of metric and scale to create a n x n distance matrix that is analyzed as a graph using a novel algorithm. Details of the algorithm
are provided in the paper cited below, by D. Baron and B. Ménard.


```julia
sequence(A, metrics=ALL_METRICS, scales=nothing)
```

or equivalently:

```julia
sequence(A)
```

as `metrics=ALL_METRICS` and `scales=nothing` are defaults. If you want to specify only one metric, you must wrap it in a 1-tuple.
e.g. to use only KL Divergence, write:

```julia
sequence(A, metrics=(KLD,), scales=nothing
```

Use the `scales` keyword to specify the number of "chunks" into which the data should be divided, as a tuple, e.g. 

```julia
sequence(A, scales=(1,3,5))
```
In the `autoscale` mode (enabled by default as `scales=nothing`) SequencerJ will find its own "best scale" based on 
running the Sequencer against a sample of columns (10% for now) and picking the scale that results in the greatest elongation.

A 1-D grid may be provided. The grid - whose deltas figure in the distance calculations - must be non-negative 
real numbers (Float16, Float32, or Float64). The grid length must equal the number of rows in A.

```julia
julia> sequence(A; metrics=(WASS1D,L2), grid=collect(0.5:0.5:size(A,1))) # grid must equal the size of A along dim 1
```

The paper that describes the Sequencer algorithm and its applications can be found 
on Arxiv: [https://arxiv.org/abs/2006.13948].
```bibtex

@misc{baron2020extracting,
    title={Extracting the main trend in a dataset: the Sequencer algorithm},
    author={Dalya Baron and Brice Ménard},
    year={2020},
    eprint={2006.13948},
    archivePrefix={arXiv},
    primaryClass={cs.LG},
    year=2020
}
"""
function sequence(A::VecOrMat{T}; 
    scales=nothing, #auto-scale by default
    metrics=ALL_METRICS,
    grid=nothing,
    weightrows=false,
    rowfn=i2weight,
    ) where {T <: Real}

    @assert all(sign.(A) .>= 0) && all(!any(isnan.(A)) && !any(isinf.(A))) "Input data cannot contain negative, NaN or infinite values."

    @info "Sequencing data with
    shape: $(size(A))
    metric(s): $(metrics)
    scale(s): $(scales)"

# We need to ensure that A is, or becomes, oriented so that
# observation data run in columns. Each column is a separate data series.
# Each row typically represents a distinct point in time or space where the
# observed data were captured.
    # take a vector of vectors and cat it into a matrix
    A = A isa Vector ? hcat(A...) : A

    # package the arguments to be consumable as iterators
    tuplify(x::PreMetric) = tuple(x)
    tuplify(x) = tuple(x...)
    
    # replace zeros and Infs with values that work for the math
    clamp!(A, eps(), typemax(eltype(A)))
    @debug "before ensure" grid
    grid = ensuregrid!(A, grid)
    @debug "after ensure" grid

    metrics = (metrics !== nothing) ? tuplify(metrics) : ALL_METRICS
    if scales === nothing
        scales = _bestscale(A, metrics, grid, (weightrows, rowfn))
        @info "After autoscaling, best scale = $(scales)..."
    end
    scales = tuplify(scales)

    return _sequence(A, scales, metrics, grid, weightrows, rowfn)
end

"The scales used in the autoscale option"
const FIB = [1,2,3,5,8,13,21,34,55,89,144,233,377]

"Sample columns of A and run the Sequencer algorithm against this subspace to identify the 
scale at which elongation is greatest. Scales are members of the Fibonacci series, e.g. 1,2,3,5,8,13..."
function _bestscale(A, m, g, params)
    M,N = size(A)
    Ŋ = N ÷ 10  # nb of column samples. CAN THIS BE SMALLER?
    sampidx = sample(collect(axes(A,2)), (Ŋ,), ordered=true)
    Asamp = view(A, :, sampidx)
    minscale = 1
    maxscale = M ÷ 10 # no fewer than 10 rows = observations
    testscales = filter(i->i < maxscale, FIB)
    bestscale = 1
    max_η = 0
    for s in testscales
        r = _sequence(Asamp, s, m, g, params...)
        η = elong(r)
        if η > max_η
            max_η = η
            bestscale = s
        end
    end
    return bestscale
end

"""
The full algorithm. For internal use only.
"""
function _sequence(A, scales, metrics, grid, weightrows, rowfn)

    # rows M and columns N in A
    M, N = size(A)

    MST_all = []
    η_all = []
    EOSeg = Dict{Tuple,Any}() # list of elongation and orderings (BFS,DFS), one per Segment    
    EOAlgScale = Dict{Tuple,Any}() # elongation and orderings for the cumulated weighted distances
    Wr = ones(M) # identity
    rowseq = collect(1:M) # row indices for applying row weight to column vectors
    if weightrows
        # force grid back to nothing here so that row sequencing uses its own grid
        # use all metrics at scale 1 for rows
        r = sequence(permutedims(A), scales=(1,), metrics=ALL_METRICS, grid=nothing, weightrows=false);
        # get the optimal ordering for rows
        # method call seems to work only when fully qualified.
        rowseq = SequencerJ.order(r)
        # weights are simply the reciprocal. Maybe look at a different formula? Linear? Wr = (1 .- rowseq ./ M)
        Wr = rowfn(rowseq)
    end
    #row weight index; use this below to reorder Wr before chunking
    rwidx = sortperm(rowseq)

    @inbounds for k in metrics
        alg = nothing
        for l in scales

            # summary distance matrix for segments
            Dklms = zeros(N, N)
            # Elongation per chunk
            ηs = []
            # BFS sequence per chunk
            orderings = []
            # split the data A, grid and row weights Wr into chunks
            S, G, W = _splitnorm(A, grid, Wr[rwidx], l)

            # Each m row in S contains n segments of data,
            # one for each of n data series
            tt = @elapsed for i in eachindex(S, G, W)
                # m is a chunk of the input data matrix
                # local grid with the same dim 1 dimension as m
                m, lgrid = S[i], G[i]
                @debug "after splitting" size(m), size(lgrid)
                # Some algorithms can be performed on arbitrary real-valued grids
                @debug "k" k
                if k in (EMD, Energy)
                    @debug "adding local grid to alg" k lgrid
                    alg = k((lgrid, lgrid))
                else
                    alg = k
                end
                @debug "alg" alg
                # Weight the columns by the appropriate row weights (default = 1)
                # @debug "W[i] m[:,1:end]" W[i] m[:,1:end]
                mp = W[i] .* m[:,1:end]
                # @debug "W[i] .* m[:,1:end]" mp
                m .= mp
                # map!( (c) -> W[i] .* c, m, collect(eachcol(m)))

                # All the heavy lifting happens here, in the distance calculations
                Dklm = abs.(pairwise(alg, m; dims=2)) .+ ϵ

                # Measure our per-metric, per-scale, per-segment distance matrix
                MSTklm, _, η, bfso = _measure_dm(Dklm)
                # weight the distance matrix by its elongation factor
                Dklm_e = η .* Dklm
                # ensure we stay on the playing field...
                clamp!(Dklm_e, eps(), typemax(eltype(Dklm_e)))
                # update the summary matrix
                Dklms .+= Dklm_e
                # save it for later...don't run away and let me down
                push!(ηs, η)
                push!(orderings, bfso)
            end

            if length(Dklms) < 1
                @warn("Unable to create distance matrices from the given data.")
                return nothing
            end

            # Perform the same MST-based analysis on the weighted results over all segments
            Dkl = Dklms / sum(ηs)

            # measure the per-metric, per-scale average distance matrix
            # minimum spanning tree, elongation, and the breadth-first search tree
            MSTkl, _, ηkl, BFSkl = _measure_dm(Dkl)

            push!(MST_all, MSTkl) # All MSTs (1 per alg,s combo)
            push!(η_all, ηkl) # All elongations (1 per alg,s combo)

            # Store these as intermediate results
            EOSeg[(alg, l)] = (ηs, orderings)
            EOAlgScale[(alg, l)] = (ηkl, BFSkl)

            @info "$(k) at scale $(l): η = $(@sprintf("%.4g", ηkl)) ($(@sprintf("%.2g", tt))s)"        
        end
    end
    
    # Create a sparse N x N proximity matrix to be filled with MST elongation-weighted edge
    # distances
    Pw = _weighted_p_matrix(N, MST_all, η_all)
    # Invert the proximity matrix to get a final distance matrix
    # w/ zeros on the diagonal
    D = inv.(Pw)
    # final minimum spanning tree for analysis
    mstD, stidx, ηD, bfstD = _measure_dm(D)
    # convert the final sequence from tree form into an ordered vector of vertices.
    order = unroll(bfstD, stidx)
    @info "Final: η = $(@sprintf("%.4g", ηD)) sequence = $(prettyp(order,5))"
    return SequencerResult(EOSeg, EOAlgScale, D, mstD, ηD, order)
end

"""
Evaluate a given distance matrix and return its MST graph, elongation, and optimal ordering of vertices. Each vertex
in the result corresponds to a column of data in D.

Internal use only.
"""
function _measure_dm(D::AbstractMatrix)
    g = _mst(D) # minimum spanning tree for D
    stidx = _leastcentralpt(g) # Least central point of MST
    η = elongation(g, stidx)
    bfst = bfs_tree(g, stidx)
    return g, stidx, η, bfst
end

"""
Ensure that the size of the 1-D grid, if provided, is compatible with the data in `A`.
Create a grid if one was not provided.
"""
function ensuregrid!(A, grid=nothing)::AbstractVector
    if isnothing(grid)
        grid = float(collect(axes(A, 1)))
    elseif !(eltype(grid) isa AbstractFloat)
        grid = float(grid)
    end
    @assert length(grid) == size(A, 1) "Grid length must match dims 1 (row count) of A. Got grid:$(length(grid)) and A:$(size(A, 1))"
    grid
end

"Create a N x N proximity matrix from weighted edges of the given graph(s).
Edge weights are provided as η coefficient(s) >= 1.
For internal use only.
"
function _weighted_p_matrix(N, graphs, ηs)
    P = spzeros(N, N)
# fill the proximity matrix with elongation-weighted reciprocal edge weights, treated as distances.
    @inbounds for idx in eachindex(graphs, ηs)
        g = graphs[idx]
        W = LightGraphs.weights(g)
        η = ηs[idx]
        for e in edges(g)
            i, j = src(e), dst(e)
            d = W[i,j] # weight as distance
            # proximity is the inverse of distance
            P[i,j] = P[j,i] += η * 1.0 / d
        end
    end
    # All P[i] were elongated using an η. Now divide by Ση to get the elongation-weighted average.
    Pw = P / sum(ηs)
    # Put Inf on the 1,1 diagonal to form a true proximity matrix
    Pw .+= sparse(Matrix(Inf * I, size(Pw)...))
    return Pw
end

"Return the index of the least central node in the given graph."
_leastcentralpt(g) = last(findmin(closeness_centrality(g)))

"Return the half-length or mean of the given path lengths (distances)."
_halflen(paths) = mean(paths)

"Return the half-width or mean count of unique values in the given path lengths."
_halfwidth(paths) = mean(count(v -> v == k, paths) for k in unique(paths)) / 2

"""`julia`

`elongation(g, startidx)`

Returns the ratio of the graph g's half-length (mean of path distances) over the 
half-width, defined as the mean count of shortest paths from the center node of a
minimum spanning tree over the graph. This function calls the `dijkstra_shortest_paths` 
function in `LightGraphs`.
"""
function elongation(G, startidx)
    spaths = dijkstra_shortest_paths(G, startidx, LightGraphs.DefaultDistance())
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
# TODO Ask Dalya about this
    g = SimpleWeightedGraph(issymmetric(dm) ? dm : Symmetric(dm))
# Minimum Spanning Tree using the Kruskal algorithm
    mst = kruskal_mst(g; minimize=true)
    Ii = src.(mst)
    Ji = dst.(mst)
    Vi = weight.(mst)
# Populate a result distance matrix with the MST nodes and edge weights and convert it to a graph
    S = sparse(Ii, Ji, Vi, size(dm)...)
    return issymmetric(S) ? SimpleWeightedGraph(S) : SimpleWeightedGraph(Symmetric(S))
end


"Split the given `m x n` matrix and grid into `scale` parts along the column dimension such
that the resulting matrices (chunks) are approximately of dimension `k x n` where `k ≊ m / scale`"
function _splitnorm(A::AbstractMatrix{T}, grid, rw, scale) where {T <: Real}

    @assert scale <= length(A[:,1]) "Scale ($(scale)) cannot be larger than the number of data elements ($(length(A[:,1])))."

    M = size(A, 1)
    # break up each column of A into  chunks
    chunklen = cld(M, scale)
    slices = []
    grids = []
    weights = []

    chunkstrt = 1:chunklen:M # the starting indices for each chunk

    for i in chunkstrt
        ii = i + chunklen - one(i) # the end index for the chunk
        ii = clamp(ii, one(ii), M)  # keepin' it real....
        S = A[i:ii, :] # a subset of rows, all columns
        G = grid[i:ii]
        W = rw[i:ii]
        # now normalize each column of the slice
        Σslice = sum(S, dims=1)
        push!(slices, S ./ Σslice)
        push!(grids, G)
        push!(weights, W)
    end
    return slices, grids, weights
    end

EOSeg = Dict{Tuple,Any}() # Dictionary of per-segment intermediate elongation and sequence results
EOAlgScale = Dict{Tuple,Any}() # elongation and orderings for the cumulated weighted distances
