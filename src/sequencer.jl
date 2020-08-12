

const L2 = Distances.SqEuclidean()
const WASS1D = EMD()
const KLD = Distances.KLDivergence()
const ENERGY = Energy()
const ALL_METRICS = (L2, WASS1D, KLD, ENERGY)

# Smallest weight/distance (instead of 0 or Inf)
const ϵ = 1e-6

# dictionaries to hold intermediate results, keyed by (metric, scale)
# list of elongation and orderings (BFS,DFS), one per Segment
d_e_o = Dict()
# elongation and orderings for the cumulated weighted distances
d_w = Dict()

"""
# variables from the python code
# INPUT
# self.grid = grid
# self.objects_list = objects_list
# self.estimator_list = estimator_list
# self.scale_list

# OUTPUT
# self.weighted_elongation_and_sequence_dictionary = None
# self.final_mst_elongation = None
# self.final_mst = None
# self.final_sequence = None
"""

"""
    sequence(A, scales=(1,4), metrics=(:all))

documentation
"""
function sequence(A::VecOrMat{T}; scales=(1,4), metrics=ALL_METRICS, grid=nothing) where T

    combos = [(m,s) for m in metrics, s in scales]

# We need to ensure that A is, or becomes, oriented so that
# observation data run in columns. Each column is a separate data series.
# Each row typically represents a distinct point in time or space where the
# observed data were captured.
#
    @assert all(sign.(A) .>= 0) && all(!any(isnan.(A)) && !any(isinf.(A))) "Input data cannot contain NaN or infinite values."

    # take a vector of vectors and cat it into a matrix
    A = A isa Vector ? hcat(A...) : A
    @debug "A " A

    @debug "scales,metrics,grid" scales metrics grid
    # replace zeros with epsilon (itsy bitsy teeny tiny value)
# DONE See if this map step can be removed safely. Answer: No, it must remain in place!
    map!(v->v ≈ 0. ? v+eps() : v, A, A)

    # create a sensible grid if one was not provided
    if isnothing(grid)
        grid = float(collect(axes(A,1)))
    elseif !(eltype(grid) isa AbstractFloat)
        grid = float(grid)
    end
    @assert length(grid) == size(A,1) # down the column...
    @debug "grid after creation/floatation" grid

    MST_all = []
    η_all = Float64[]
    # Dkls = []
    for (alg,s) in combos
    # build a master dictionary of distance matrices for all metrics and scales
        # Dijk = Dict{Tuple,Vector{Array{T}}}()
        @info "Metric $(alg) at scale $(s)..."

        S, G = _splitnorm(s, grid, A)
        @debug "S after split" S

        Dklms = zeros(size(A,2), size(A,2))
        ηs = [] # Elongations per chunk
        orderings = [] # BFS, DFS orderings per chunk

        # Each m row in S contains n segments of data,
        # one for each of n data series
        @inbounds for i in eachindex(S, G) # SEGMENTS
            #convert to a matrix
            m = S[i] # m is a matrix already
            localgrid = G[i] # a vector
            @debug "matrix for distance calcs" m
            # m = hcat(r...)'
            if alg in (EMD, Energy)
                alg = alg((localgrid, localgrid))
            end
            Dklm = abs.(pairwise(alg, m; dims = 2)) .+ ϵ
            @debug "Dklm" Dklm

            MSTklm = _mst(Dklm) #MST
            startidx = _startindex(MSTklm)

            η = elongation(MSTklm, startidx)
            η = isnan(η) || isinf(η) ? typemin(η) : η
            @debug "elongation for Dklm" η
            push!(ηs, η)

# weight the distance matrix by its elongation factor
            Dklm_e = η .* Dklm
            clamp!(Dklm_e, eps(), typemax(eltype(Dklm_e)))
            # map!(v -> isinf(v) ? typemax(v) : v, , Dklm_e)
            @debug "Dklm_e after elongation" Dklm_e
            Dklms .+= Dklm_e

            bfso, dfso = _b_d_order(MSTklm, startidx)
            @debug "BFS,DFS" bfso dfso
            push!(orderings, (bfso,dfso))
        end

        if length(Dklms) < 1
            @warn("Unable to create distance matrices from the given data.")
            return nothing
        end

        @debug "Dklms after accumulation" Dklms
        Dkl = Dklms / sum(ηs) # Weighted average over all chunks for <alg,s>
        @debug "elongation-weighted Dkl" Dkl
        MSTkl = _mst(Dkl)
        @debug ""
        push!(MST_all, MSTkl) # All MSTs (1 per alg,s combo)
        stidx = _startindex(MSTkl) # Least central point of averaged MST
        ηkl = elongation(MSTkl, stidx)
        @debug "elongation of weighted Dkls" ηkl
        push!(η_all, ηkl) # All elongations (1 per alg,s combo)

        BFSkl = bfs_tree(MSTkl, stidx)

        # Store these as intermediate results
        global d_e_o[(alg,s)] = (ηs, orderings)
        global d_w[(alg,s)] = (ηkl, BFSkl)
    end

    # Weighted average of all metrics, scales and chunks (klm)
    # Dall = sum(Dkls, dims=3) / sum(η_all)
    # @debug Dall
    # Dall seems not to be used by the rest of the algorithm!!

    # N = maximum(nv.(MST_all)) # vector of Graphs
    # @debug "maximum of nv.(MST_all)" N
    N = size(A,2)

    # A sparse proximity matrix to be filled with MST elongation-weighted edge
    # distances (which are actually weights as well)
    P = spzeros(N,N)

    @inbounds for idx in eachindex(MST_all, η_all)
        g = MST_all[idx]
        W = LightGraphs.weights(g)
        η = η_all[idx]
        for e in edges(g)
            i,j = src(e),dst(e)
            d = W[i,j]
            P[i,j] = P[j,i] += η * 1.0 / d
        end
    end
    # All P[i] were elongated using an η. Now divide by Ση to get the elongation-weighted average.
    Ση = sum(η_all)
    Pw = P / Ση
    # Put Inf on the 1,1 diagonal to form a true proximity matrix
    Pw .= sparse(Matrix(Inf*I, size(Pw)...))
    @debug "Amazing! Pw Pw Pw!" Pw

    # Invert the proximity matrix to get a final distance matrix, w/ zeros
    # on the diagonal
    D = similar(Pw)
    D .= 1 ./ Pw
    @debug "D" D

    # final minimum spanning tree for analysis
    MSTD = _mst(D) #MST
    @debug "final MSTD" MSTD
    stidx = _startindex(MSTD) # Least central point of averaged MST
    # Here's that elongation thing again...
    ηD = elongation(MSTD, stidx)
    @debug "elongation of MSTD " ηD

    BFSD = bfs_tree(MSTD, stidx)
    @debug "final BFSD" BFSD

    return MSTD, ηD, BFSD
end

"""
Returns the index of the least central node in the given graph.
"""
function _startindex(g)
    @debug "starting cc calculation for graph $(g)..."
    cc = closeness_centrality(g) # a vector of measures, one per node in the graph
    @debug "closeness_centrality for g $(g)" cc
    _, minidx = findmin(cc) # value and index of minimum val in C
    return minidx
end


halflen(paths) = mean(paths)

halfwidth(paths) = mean(count(v->v==k, paths) for k in unique(paths)) / 2

"""`julia`

`elongation(g, startidx)`

Returns the ratio of the graph half-length (mean of path distances) over the 
half-width, defined as the mean count of shortest paths from the center node of a
minimum spanning tree over the graph.
"""
function elongation(g, startidx)
    spaths = dijkstra_shortest_paths(g, startidx)
    @debug "spaths.dists" spaths.dists
    # remove those annoying zeros and Infs that ruin it for everyone else
    hlen = clamp(halflen(spaths.dists), eps(), floatmax(Float32))
    hwid = clamp(halfwidth(spaths.dists), eps(), floatmax(Float32))
    @debug "halflen, halfwidth" hlen hwid
    return hlen / hwid + 1
end

# BREADTH FIRST: Uses alternate outneighbors implementation; see below
_b_d_order(g, minidx) = (bfs_tree(g, minidx), dfs_tree(g, minidx))

"""
Returns a weighted graph of the MST of the given distance matrix,
calculated using the Kruskal algorithm.
"""
function _mst(dm::AbstractMatrix)
    # Construct a graph from the distance matrix.
    # By forcing an unsymmetric distance matrix into a Symmetric one
    # we are losing 1/2 (the lower half) of the distances
    # Should we construct a separate graph for the lower half as well?

    map!(v -> v ≈ 0. ? ϵ : v, dm, dm)
    g = issymmetric(dm) ? SimpleWeightedGraph(dm) : SimpleWeightedGraph(Symmetric(dm))
    @debug "graph from distance matrix ne nv" g ne(g) nv(g)
    @debug "minimum weight value" minimum(LightGraphs.weights(g))

# Minimum Spanning Tree using the Kruskal algorithm
    mst = kruskal_mst(g; minimize=false) # calculate the maximum tree instead of the minimized version
    @debug "kruskal MST" mst
    Ii = src.(mst)
    Ji = dst.(mst)
    Vi = weight.(mst)
    @debug "I J V" Ii Ji Vi

    S = sparse(Ii, Ji, Vi, size(dm)...)
    @debug "big beautiful S symmetric?" S issymmetric(S)

# It really needs to be possible to construct a SWG using an Edge iterator...
    # mstg = SimpleWeightedGraph(length(mst))
    # for e in mst
    #     add_edge!(mstg, e)
    # end
    # @debug "weighted graph from DM ne nv" mstg ne(mstg) nv(mstg) 
    # mstg
    g = issymmetric(S) ? SimpleWeightedGraph(S) : SimpleWeightedGraph(Symmetric(S))
    return SimpleWeightedGraph(g)
end


"""
Calculates a distance matrix for A for the given
algorithm. The scale and returns a Dictionary of DMs.
"""
function _splitnorm(scale, grid, A)
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


import LightGraphs: bfs_parents, _bfs_parents

bfs_parents(g::AbstractGraph, s::Integer; dir = :out) =
    (dir == :out) ?
    LightGraphs._bfs_parents(g, s, outneighbors_ranked) :
    LightGraphs._bfs_parents(g, s, inneighbors)

"""
default order = :desc (highest weighted edges first). Any other value
means order ascending.
"""
function outneighbors_ranked(g, v; order=:asc)
    alln = collect(outneighbors(g,v))
    # @show "Before $(alln)"
    w = LightGraphs.weights(g)
    T = eltype(w)
    rw = T[]
    for n in alln
        push!(rw, w[v,n])
    end
    idx = zeros(Int,length(rw))
    sortperm!(idx, rw; rev = (order == :desc))
    return alln[idx]
end



########################################################################################################
        ######### STEP 2: order the spectra based on the different distance matrices, and measure    ###########
        #########         weights using the MST elongations.                                         ###########
        #########         Produce weighted distance matrix per scale and estimator.                  ###########
        ########################################################################################################


        ########################################################################################################
        ######### STEP 3: use the elongations of the weighted distance matrices sequences            ###########
        #########         to build proximity matrices, then convert them to distance matrices,       ###########
        #########         and obtain the final BFS and DFS sequences.                                ###########
        ########################################################################################################


        ########################################################################################################
        ######### STEP 4: save the final BFS and DFS sequences, their final elongation, and          ###########
        #########         the sparse distance matrix that was used to obtain these.                  ###########
        ########################################################################################################
