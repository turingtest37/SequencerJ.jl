module SequencerJulia

using Reexport
@reexport using Distances
using StatsBase
using SparseArrays
using LightGraphs
using SimpleWeightedGraphs
using Logging
using LinearAlgebra


include("distancemetrics.jl")
include("sequencer.jl")

# functions
export emd, energy, cdf_distance, sequence

# types
export EMD, Energy

export L2, WASS1D, KLD, ENERGY, ALL_METRICS

# temp
export _splitnorm, _mst, _startindex, elongation, _b_d_order

end
