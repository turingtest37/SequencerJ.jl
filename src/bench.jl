using SequencerJulia
using BenchmarkTools
using Images


imgmed = load(joinpath(@__DIR__,"..","resources","bread.jpeg"));
AM = convert(Matrix{Float32}, float32.(Gray.(imgmed)));
println("Matrix size = $(size(AM))")

for alg in ALL_METRICS
    println(alg)
    @btime abs.(pairwise($alg, $AM; dims = 2)) .+ eps()
end
