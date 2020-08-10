using SequencerJulia
using BenchmarkTools
using Images


imgmed = load(joinpath(@__DIR__,"..","resources","bread.jpeg"));
AM = convert(Matrix{Float32}, float32.(Gray.(imgmed)));
println("Matrix size = $(size(AM))")

AM ./= sum(AM, dims=1)

# AMR = permutedims(AM,[2,1])

# AM = _splitnorm(1, collect(axes(AM,1)), AM)
# Profile.@profile begin
for alg in ALL_METRICS
    println(alg)
    @btime abs.(pairwise($alg, $AM; dims = 2)) .+ eps()
end
# end

# Profile.print(; noisefloor=2.0)
#
# open("/tmp/profile.txt", "w") do s
#     Profile.print(IOContext(s, :displaysize => (24, 500)))
# end
