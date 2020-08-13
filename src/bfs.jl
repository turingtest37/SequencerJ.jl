"""
Overload LightGraphs.bfs_parents to provide an alternative outneighbors function.

See [`LightGraphs.bfs_parents`](@Ref)
"""
bfs_parents(g::AbstractGraph, s::Integer; dir = :out) =
    (dir == :out) ?
    LightGraphs._bfs_parents(g, s, outneighbors_ranked) :
    LightGraphs._bfs_parents(g, s, inneighbors)

"""
For the given graph and vertex, return the list of the vertex's out neighbors 
ranked by edge weight (ascending by default).
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
