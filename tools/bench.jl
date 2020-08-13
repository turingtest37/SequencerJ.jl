using SequencerJulia
using BenchmarkTools
using Images


runbench(A, alg) = @benchmark abs.(pairwise($alg, $A; dims = 2)) .+ eps()

imgmed = load(joinpath(@__DIR__,"..","resources","bread.jpeg"));
AM = convert(Matrix{Float32}, float32.(Gray.(imgmed)));
println("Matrix size = $(size(AM))")

AM ./= sum(AM, dims=1)

for alg in ALL_METRICS
    runbench(AM, alg)
end


# AMR = permutedims(AM,[2,1])

# AM = _splitnorm(1, collect(axes(AM,1)), AM)
# Profile.@profile begin
# end

# Profile.print(; noisefloor=2.0)
#
# open("/tmp/profile.txt", "w") do s
#     Profile.print(IOContext(s, :displaysize => (24, 500)))
# end
