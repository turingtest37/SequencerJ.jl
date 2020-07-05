module SequencerJulia

using Distances
using StatsBase


# functions
export l2, kl, wasserstein1d, emd, energy, cdf_distance

# types
export EMD

include("sequencer.jl")
include("distancemetrics.jl")

end
