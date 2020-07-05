using BenchmarkTools
using SequencerJulia

"""
``
julia> objs = [[1,1,1], [2,3,4], [5,6,7], [8,9,10]]
4-element Array{Array{Int64,1},1}:
 [1, 1, 1]
 [2, 3, 4]
 [5, 6, 7]
 [8, 9, 10]

julia> datalen = maximum(first.(size.(objs)))
3

julia> A = reshape(collect(Iterators.flatten(objs)), datalen, :)
3×4 Array{Int64,2}:
 1  2  5   8
 1  3  6   9
 1  4  7  10

julia> permutedims(A)
4×3 Array{Int64,2}:
 1  1   1
 2  3   4
 5  6   7
 8  9  10

``
"""
function method_a(A)
    datalen = maximum(first.(size.(A)))
    return permutedims(reshape(collect(Iterators.flatten(A)), datalen, :))
end


method_b(A) = hcat(A...)'

method_c(A) = permutedims(hcat(A...))

A = [[1,1,1], [2,3,4], [5,6,7], [8,9,10]]

B = [rand(UInt16,1000) for _ in 1:1000];

# for f in (:method_a, :method_b, :method_c)
#     @btime f($(A))
# end
for m in (A,B)
    for x in (:a,:b,:c)
        println("Method $(x) with array size $(size(m))")
        @eval @btime :(Symbol("method_",$x)($$m))
    end
end
