using SequencerJulia
using LightGraphs
using SimpleWeightedGraphs
using LinearAlgebra
using Statistics: mean

using BenchmarkTools
using Images

using UnbalancedOptimalTransport


a = rand(100)
Σa = sum(a)
wa = a / Σa
A = Float32.(collect(1:100))

b = rand(100)
Σb = sum(b)
wb = b / Σb
B = Float32.(collect(1:100))

a = DiscreteMeasure(wa, A)
b = DiscreteMeasure(wb, B)

cost = (x,y) -> abs(x - y)
ϵ = 0.01
D = UnbalancedOptimalTransport.KL(1)
SD = sinkhorn_divergence!(D, cost, a, b, ϵ)
π = sum(optimal_coupling!(D, cost, a, b, ϵ))

k = emd(wa, wb)

imgmed = load(joinpath(@__DIR__,"..","resources","bread.jpeg"));
AM = convert(Matrix{Float32}, float32.(Gray.(imgsmall)));

imgsmall = load(joinpath(@__DIR__,"..","resources","colony.png"));
A = convert(Matrix{Float32}, float32.(Gray.(imgsmall)));

A = A isa Vector ? hcat(A...) : A
# if isnothing(grid)

map!(v->v≈0. ? v+eps() : v, A, A)


mygrid = Float32.(collect(axes(A,1)))

m = alg = KLDivergence()
s = 2

len = cld(size(A,1), s)
N = cld(size(A,1), len) # Number of chunks
@debug "N (number of bins)" N

slices=[]
# range of chunk starting indices
r = 1:len:size(A,1)
for i in r
    ii = i + len - 1 # thank you, 1-based indexing
    # use the "end" keyword on the last slice to mop up the remainder
    ii = ii > size(A,1) ? :end : ii

    S = @eval A[$i:$ii, :]
    # @show size(S)
    # now normalize each column of the slice
    Σslice = sum(S, dims=1)
    # @show Σslice

    Snorm = S ./ Σslice
    # @show sum(norms, dims=1)
    push!(slices, Snorm)
end
slices
# chunks, S = _splitnorm(s, mygrid, A)

# SR = reshape(chunks, 3, :)
# SR[1,1]
# hcat(chunks...)

Dklms = Matrix{AbstractFloat}[] #maybe thi needs to be zeros
ηs = [] # Elongations per chunk
orderings = [] # BFS, DFS orderings per chunk

# r = collect(Iterators.take(eachrow(slices), 1))

# size(r[1])


#[3,:] # SEGMENTS
#convert to a matrix
# m = hcat(r...)'

m = first(slices)

Dklm = abs.(pairwise(alg, m; dims = 1)) .+ eps()
 #DISTANCE MATRIX
@debug "Dklm" Dklm

MSTklm = _mst(Dklm) #MST
startidx = _startindex(MSTklm)
@debug "start idx for MSTklm" startidx

spaths = dijkstra_shortest_paths(MSTklm, startidx)
η = elongation(MSTklm, startidx) + 1
@debug "elongation for Dklm" η
if isnan(η) || isinf(η)

end
push!(ηs, η)

# weight the distance matrix by its elongation factor
Dklm_e = η .* Dklm
@debug "Dklm_e after elongation" Dklm_e
push!(Dklms, Dklm_e)

bfso, dfso = _b_d_order(MSTklm, startidx)
@debug "BFS,DFS" bfso dfso
push!(orderings, (bfso,dfso))
# end


# dml = vec(unweighted[(m,s)])

# dml = first(collect(values(unweighted)))
# for (m,s) in combos
#     dml = vec(unweighted[(m,s)])
    # wdm, ocl, ecl = weight(dml)
    # weighted_distance_matrix, ordering_per_chunk_list, elongation_per_chunk_list = self._return_weighted_distance_matrix_for_single_estimator_and_scale(distance_matrix_list, to_return_elongation_list=True)

# end

# dm = first(collect(values(dml)))

g = issymmetric(dm) ? SimpleWeightedGraph(dm) : SimpleWeightedGraph(Symmetric(dm))

mst = kruskal_mst(g) # vector of Edges

SimpleWeightedGraph(mst)

mstg = SimpleWeightedGraph(length(mst))

for e in mst
    add_edge!(mstg,e)
end

C = closeness_centrality(mstg) # a vector of measures, one per node in the graph

(minx,minidx) = findmin(C) # value and index of minimum val in C

bfs_order = bfs_tree(g, minidx)

dfsorder = dfs_tree(mstg, minidx)

spaths = dijkstra_shortest_paths(mstg, minidx)

sdists = spaths.dists

half_len = mean(sdists)

d = sum(count(v->v==k, sdists) for k in unique(sdists))

half_width = mean(d) / 2

mst_elongation = half_len / half_width + 1

valtype()

## Or, better yet, import and overload a function in LightGraphs directly
bfsorder = bfs_tree(mstg, minidx)


collect(edges(C))
C
