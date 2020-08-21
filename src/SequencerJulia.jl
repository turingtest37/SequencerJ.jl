module SequencerJulia

using Logging
using Printf: @sprintf
using Reexport
@reexport using Distances
using SparseArrays: sparse, spzeros
using LightGraphs
using SimpleWeightedGraphs
using StatsBase: mean, ecdf, ECDF
using LinearAlgebra: issymmetric, Symmetric, I
import LightGraphs: bfs_parents, _bfs_parents
import Base: show


include("bfs.jl")
include("utils.jl")
include("distancemetrics.jl")
include("sequencer.jl")

# functions
export emd, energy, cdf_distance, sequence

# types
export EMD, Energy

export L2, WASS1D, KLD, ENERGY, ALL_METRICS

export elong_start_index, elongation, D, mst, elong, stidx, order, rollup, prettyp


end
