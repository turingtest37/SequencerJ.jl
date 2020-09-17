# API

```@meta
CurrentModule = SequencerJ
```

## Types

```@docs
SequencerResult
```

Call accessor functions using the SequencerResult object. For example, to get the final column ordering of the data, call [`order`](@ref) :
```@example
using SequencerJ #hide
A = rand(100,100);
r = sequence(A;);
order(r)
```

## Functions
```@docs

sequence
elong
order
mst
D
η

EMD()
EMD(u::AbstractVector{T}, v::AbstractVector{T}) where {T <: Real}
EMD(u,v,uw,vw)
emd

Energy()
Energy(::AbstractVector{T}, ::AbstractVector{T}) where {T <: Real}
Energy(u,v,uw,vw)
energy

elongation
ensuregrid!
unroll
prettyp
```

## Metrics

  * [`KLD`](@ref)
  * [`L2`](@ref)
  * [`ENERGY`](@ref)
  * [`WASS1D`](@ref)
  * [`ALL_METRICS`](@ref)
  * [`Euclidean()`](@ref)
  * [`KLDivergence()`](@ref)
