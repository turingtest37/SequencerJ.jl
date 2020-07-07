module SequencerJulia

using Distances
using StatsBase


# functions
export l2, kl, kldm, wasserstein1d, emd, energy, cdf_distance

# types
export EMD, Energy

import Distances: evaluate

include("sequencer.jl")
include("distancemetrics.jl")

end
