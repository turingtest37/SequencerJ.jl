"""
Overload of LightGraphs.bfs_parents that provides an alternative outneighbors function.

See [`LightGraphs.bfs_parents`](@ref)
"""
LightGraphs.bfs_parents(g::AbstractGraph, s::Integer; dir = :out) =
    (dir == :out) ?
    LightGraphs._bfs_parents(g, s, outneighbors_ranked) :
    LightGraphs._bfs_parents(g, s, inneighbors)

"""
For the given graph and vertex, return the list of the vertex's out neighbors 
ranked by edge weight between vertex and the neighbor (ascending order by default).
This function is invoked by [`bfs_parents`](@ref) and is not intended to be called directly.
    @ TODO #17 Benchmark this to see if using views would be faster.
"""
function outneighbors_ranked(g, v; order=:asc)
    alln = collect(outneighbors(g,v))
    length(alln) == 1 && return alln
    W = LightGraphs.weights(g)
    T = eltype(W)
    rw = T[]
    for n in alln
        w = W[v,n]
        push!(rw, w)
    end
    idx = Array{Int}(undef, length(rw))
    sortperm!(idx, rw; rev = (order == :desc))
    return alln[idx]
end
