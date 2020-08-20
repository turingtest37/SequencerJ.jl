"""
Print the sequence of outbound nodes of a graph starting from the given index.
pystyle implies reversing the final order and substracting 1 from the vertex number.
"""
function prettyp(G::AbstractGraph, strtidx::Int; pystyle=false)
    arrow = G isa LightGraphs.SimpleDiGraph ? "->" : "--"
    v = unroll(G, strtidx)
    println(join(string.(pystyle ? reverse(v .- 1) : v), arrow))
end

"""
Return the outbound vertices that were visited from the given starting index.
"""
function unroll(G::AbstractGraph, idx::Int; visited=Int[])
    idx in visited && return visited
    push!(visited, idx)
    N = outneighbors(G, idx)
    length(N) < 1 && length(visited) == nv(G) ? (return visited) : nothing
    for n in N
        unroll(G, n, visited=visited)
    end
    return visited
end
