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
Walk the given graph, starting from the given vertex,
returning the list of all outbound vertices that are visited.
"""
function unroll(G::AbstractGraph, idx::Int; visited=Int[])
    idx in visited && return visited
    push!(visited, idx)
    N = outneighbors(G, idx)
    length(N) < 1 && length(visited) == nv(G) ? (return visited) : nothing
    for n in N
        unroll(G, n, visited=visited)
    end
    return reverse(visited)
end

"""
Return the first and last `len` elements of the vector as a string.

Default is to print 3 elements at head and tail.

```julia
julia> v = collect(1:10);
julia> prettyp(v)

"1,2,3...8,9,10"
```

Example w/5 elements visible

```julia
julia> v = collect(1:10);
julia> prettyp(v, 5)

"1,2,3,4,5...6,7,8,9,10"
```

"""
function prettyp(v::AbstractVector, len::Int = 3)
    length(v) < 7 && return join(string.(v), ",")
    head = join(string.(v[1:len]),",")
    tail = join(string.(v[(end-(len-1)):end]),",")
    return join([head,tail],"...")
end