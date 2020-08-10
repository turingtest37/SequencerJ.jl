
########################################################################################################
        ######### STEP 1: load or calculate distance matrices for different estimators and scales    ###########
########################################################################################################

const L2_METRIC = Distances.SqEuclidean()
const EMD_METRIC = EMD()
const KLD_METRIC = Distances.KLDivergence()
const ENERGY_METRIC = Energy()
const ALL_METRICS = (L2_METRIC,EMD_METRIC,KLD_METRIC,ENERGY_METRIC)

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
    @assert all(!any(isnan.(A)) && !any(isinf.(A))) "Input data cannot contain NaN or infinite values."

    # take a vector of vectors and cat it into a matrix
    A = A isa Vector ? hcat(A...) : A
    @debug "A " A

    # replace zeros with epsilon itsy bitsy teeny tiny value
    map!(v->v ≈ 0. ? v+eps() : v, A, A)

    # create a sensible grid if one was not provided
    if isnothing(grid)
        grid = Float32.(collect(axes(A,1)))
    elseif !(eltype(grid) isa AbstractFloat)
        grid = Float32.(grid)
    end

    @assert length(grid) == size(A,1) # down the column...

    MST_all = []
    η_all = Float64[]
    # Dkls = []
    for (alg,s) in combos
    # build a master dictionary of distance matrices for all metrics and scales
        # Dijk = Dict{Tuple,Vector{Array{T}}}()
        @info "Metric $(alg) at scale $(s)..."

        S, G = _splitnorm(s, grid, A)
        @debug "S after split" S

        Dklms = zeros(size(A,2), size(A,2))#maybe thi needs to be zeros
        ηs = [] # Elongations per chunk
        orderings = [] # BFS, DFS orderings per chunk

        # Each m row in S contains n segments of data,
        # one for each of n data series
        @inbounds for i in eachindex(S, G) # SEGMENTS
            #convert to a matrix
            m = S[i] # r is a matrix already, yay!
            g = G[i]
            @debug "matrix for distance calcs" m
            # m = hcat(r...)'
            Dklm = abs.(pairwise(alg, m; dims = 2)) .+ eps()
            @debug "Dklm" Dklm

            MSTklm = _mst(Dklm) #MST
            startidx = _startindex(MSTklm)
            @debug "start idx for MSTklm" startidx

            η = elongation(MSTklm, startidx) + 1
            η = isnan(η) || isinf(η) ? typemin(η) : η
            @debug "elongation for Dklm" η
            push!(ηs, η)

# weight the distance matrix by its elongation factor
            Dklm_e = η .* Dklm
            map!(v -> isinf(v) ? typemin(v) : v, Dklm_e, Dklm_e)
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
        push!(MST_all, MSTkl) # All MSTs (1 per alg,s combo)
        stidx = _startindex(MSTkl) # Least central point of averaged MST
        ηkl = elongation(MSTkl, stidx) + 1
        @debug "elongation of weighted Dkls" ηkl
        push!(η_all, ηkl) # All elongations (1 per alg,s combo)

        BFSkl, DFSkl = _b_d_order(MSTkl, stidx)

        # Store this crap for later
        global d_e_o[(alg,s)] = (ηs, orderings)
        global d_w[(alg,s)] = (ηkl, (BFSkl, DFSkl))
    end

    # Weighted average of all metrics, scales and chunks (klm)
    # Dall = sum(Dkls, dims=3) / sum(η_all)
    # @debug Dall
    # Dall seems not to be used by the rest of the algorithm!!

    N = maximum(nv.(MST_all)) # vector of Graphs
    @debug "maximum of nv.(MST_all)" N
    # ???? Why ??
    Ση = sum(η_all) # vector of Floats

    P = zeros(typeof(Ση),N,N)
    @debug "P initial" P

    @inbounds for idx in eachindex(MST_all, η_all)
        g = MST_all[idx]
        W = LightGraphs.weights(g)
        η = η_all[idx]
        for e in edges(g)
            i,j = src(e),dst(e)
            d = W[i,j]
            P[i,j] = P[j,i] += η * 1.0 / (d + eps())
        end
    end
    Pw = P / Ση
    di = [CartesianIndex(i,i) for i in 1:size(Pw,1)]
    Pw[di] .= Inf
    @debug "amazing" Pw
    ## end of Step 2

    ## Start Step 3
    D = similar(Pw)
    D .= 1 ./ Pw

    MSTD = _mst(D) #MST
    stidx = _startindex(MSTD) # Least central point of averaged MST
    ηD = elongation(MSTD, stidx) + 1
    @debug "elongation of MSTD " ηD

    BFSD, DFSD = _b_d_order(MSTD, stidx)

    return MSTD, ηD, BFSD, DFSD
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


halflen(paths) = mean(paths) + eps()

halfwidth(paths) = (mean(count(v->v==k, paths) for k in unique(paths)) + eps()) / 2

function elongation(g, startidx)
    spaths = dijkstra_shortest_paths(g, startidx)
    @debug "spaths.dists" spaths.dists
    # remove those annoying Infs that ruin it for everyone else
    # D = map(v -> clamp(v, eps(), typemax(Float32)), spaths.dists)
# elongation ratio over the shortest paths from the center node of a
# minimum spanning tree of the given distance matrix
    hlen = halflen(spaths.dists)
    hwid = halfwidth(spaths.dists)
    hlen = clamp(hlen, eps(), floatmax(Float32))
    hwid = clamp(hwid, eps(), floatmax(Float32))
    @debug "halflen, halfwidth" hlen hwid
    return hlen / hwid
end

function _b_d_order(g, minidx)
    bfsorder = bfs_tree(g, minidx) # BREADTH FIRST: Uses alternate outneighbors implementation; see below
    dfsorder = dfs_tree(g, minidx) # DEPTH FIRST
    return bfsorder, dfsorder
end

"""
Returns a weighted graph of the MST of the given distance matrix,
calculated using the Kruskal algorithm.
"""
function _mst(dm::AbstractMatrix)
    # Construct a graph from the distance matrix.
    # By forcing an unsymmetric distance matrix into a Symmetric one
    # we are losing 1/2 (the lower half) of the distances
    # Should we construct a separate graph for the lower half as well?
    g = issymmetric(dm) ? SimpleWeightedGraph(dm) : SimpleWeightedGraph(Symmetric(dm))

# Minimum Spanning Tree using the Kruskal algorithm
    mst = kruskal_mst(g) # vector of Edges
    @debug "kruskal MST" mst

# It really needs to be possible to construct a SWG using an Edge iterator...
    mstg = SimpleWeightedGraph(length(mst))
    for e in mst
        add_edge!(mstg, e)
    end
    @debug "weighted graph from DM" mstg
    mstg
end


"""
Calculates a distance matrix for A for the given
algorithm. The scale and returns a Dictionary of DMs.
"""
function _splitnorm(scale, grid, A)
    @debug "scale" scale
    # @show grid
    # @show A¡
    @assert scale <= length(A[:,1]) "Scale ($(scale)) cannot be larger than the number of data elements ($(length(A[:,1])))."

    # break up each column of A into N chunks, where N = scale
    chunklen = cld(size(A,1), scale)
    @debug "chunklen" chunklen

    N = cld(size(A,1), chunklen) # Number of chunks
    @debug "N (number of bins)" N

    slices=[]
    grids=[]
    # range of chunk starting indices
    chunkrange = 1:chunklen:size(A,1)
    for i in chunkrange
        ii = i + chunklen - 1 # thank you, 1-based indexing
        # use the "end" keyword on the last slice to mop up the remainder
        ii = clamp(ii, 1, size(A,1))
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
