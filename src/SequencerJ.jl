module SequencerJ

using Reexport
@reexport using Distances
@reexport using Logging
using Printf: @sprintf
using SparseArrays: sparse, spzeros
using LightGraphs
using SimpleWeightedGraphs
using StatsBase: mean, ecdf, ECDF, sample
using LinearAlgebra: issymmetric, Symmetric, I


import LightGraphs: bfs_parents, _bfs_parents
import Base: show


include("bfs.jl")
include("utils.jl")
include("distancemetrics.jl")
include("sequencer.jl")

# functions
export sequence, D, mst, elong, order, prettyp, loss, ensuregrid!
export emd, energy

# types
export SequencerResult, EMD, Energy

# constants
export L2, WASS1D, KLD, ENERGY, ALL_METRICS

end