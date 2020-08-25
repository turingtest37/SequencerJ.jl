
"Smallest allowed graph weight/distance (instead of 0, which does not play well with sparse graphs.)"
const ϵ = 1e-6

"""

    EMD

Earth Mover's Distance, a.k.a. 1-p Monge-Wasserstein distance

"""
struct EMD <: SemiMetric
    grids::Union{Nothing,Tuple}
end

"EMD distance with with no grid provided. Grids default to axes(u,1) and axes(v,1)."
EMD() = EMD(nothing)

"""

    EMD(u::AbstractVector{T}, v::AbstractVector{T}) where {T <: Real}

Calculate the Earth Mover Distance (EMD) a.k.a the 1-Wasserstein distance
between the two given vectors, accepting a default grid. `u` and `v` are treated as 
weights on the grid. The default grid is equal to the first axis of `u` and `v`.

This function is not intended to be called directly. Instead, use WASS1D.
"""
function EMD(u::AbstractVector{T}, v::AbstractVector{T}) where {T <: Real}
    gridu = collect(first(axes(u)))
    gridv = collect(first(axes(v)))
    EMD(gridu, gridv, u, v)
end

"""

    EMD(u,v,uw,vw)

Calculate the Earth Mover's Distance using an explicit grid. This method is not intended to be called directly.
Instead, specify a grid in the call to sequence, with the WASS1D constant.

```julia-repl

julia> A = rand(50, 100)
50×100 Array{Float64,2}:
[...]
julia> g = collect(0.5:0.5:div(size(A,1),2))
50-element Array{Float64,1}:
   0.5
   1.0
   1.5
  [...]
  24.0
  24.5
  25.0
julia> sequence(A; grid = g)
┌ Info: Sequencing data with
│     shape: (50, 100)
│     metric(s): (SqEuclidean(0.0), EMD(nothing), KLDivergence(), Energy(nothing))
└     scale(s): (1, 2, 4)
[ Info: SqEuclidean(0.0) at scale 1: η = 5.214 (3.4s)
[...]
```
"""
function EMD(u,v,uw,vw)
    ndims(u) == 1 && ndims(u) == 1 || error("Only 1-d data vectors are supported.")
    length(u) == length(uw) || error("u and u weights must have same length. Got u $(length(u)) and uw $(length(uw))")
    length(v) == length(vw) || error("v and v weights must have same length. Got v $(length(v)) and vw $(length(vw))")
    cdf_distance(ecdf(float(u); weights = float(uw)), ecdf(float(v); weights = float(vw)), 1)
end

"Same as `EMD(u,v)` but using Distances-style runnable type syntax."
function (m::EMD)(u,v=u)
    m.grids === nothing ? EMD(u,v) : EMD(m.grids...,u,v)
end

"Same as `EMD(u,v,uw,vw)` but using Distances-style runnable type syntax."
(::EMD)(u,v,uw,vw) = EMD(u,v,uw,vw)

"""
Convenience method for `EMD(u,v)`.

See [`EMD(u,v)`](@ref)
"""
emd(u::AbstractArray, v = u) = EMD(u, v)

"""
Convenience method for `EMD(u,v,uw,vw)`.

See [`EMD(u,v,uw,vw)`](@ref)
"""
emd(u,v,uw,vw) = EMD(u,v,uw,vw)


# ******** ENERGY *********

"""
    Energy

Energy distance as defined by Székely. An explicit grid may be provided. Default is 
axes(u,1) where 

[Energy distance](https://en.wikipedia.org/wiki/Energy_distance)

"""
struct Energy <: SemiMetric
    grids::Union{Nothing,Tuple}
end

"Default constructor with no specified grid."
Energy() = Energy(nothing)

"Enables dispatching on the type directly"
function (m::Energy)(u,v)
    m.grids === nothing ? Energy(u,v) : Energy(m.grids...,u,v)
end

"Direct dispatch on an object of type Energy."
(::Energy)(u,v,uw,vw) = Energy(u,v,uw,vw)


"""

    Energy(u,v,uw,vw)

Calculate Székely's energy distance between the two given vectors, `u` and `v`, 
whose weights `uw`, `vw` are treated as a empirical cumulative distribution function (CDF).
`u` and `v` must have the same length, respectively, as `uw` and `vw`. 
"""
function Energy(u,v,uw,vw)
    length(u) ==  length(uw) || error("u and u weights must have same length. Got u $(length(u)) and uw $(length(uw))")
    length(v) ==  length(vw) || error("v and v weights must have same length. Got v $(length(v)) and vw $(length(vw))")
    ndims(u) == 1 && ndims(v) == 1 || error("Only 1-d data vectors are supported.")
    sqrt(2) * cdf_distance(ecdf(float(u); weights = float(uw)), ecdf(float(v); weights = float(vw)), 2)
end


"""

    Energy(u,v)

Calculate Székely's energy distance between the two given vectors, accepting a default grid.
`u` and `v` are treated as weights on the grid. The default grid is equal to the first axis of `u` and `v`.
"""
function Energy(u,v)
    gridu = collect(first(axes(u)))
    gridv = collect(first(axes(v)))
    Energy(gridu,gridv,u,v)
end

"""
Convenience function for Energy with a supplied grid.

See [`Energy(u,v,uw,vw)`](@ref)
"""
energy(u, v, uw, vw) = Energy(u, v, uw, vw)

"""
Convenience function for Energy with a default unit grid.

    See [`Energy(u,v)`](@ref)
"""
energy(u::AbstractVector, v=u) = Energy(u,v)

"""
Returns the statistical distance between two cumulative
distributions u et v. With dimension parameter p = 1,
the distance is the 1 Wasserstein (or Earth Mover's) Distance.
With p=2, the distance is Energy.

References (from the stats.py source code)
    ----------
    [1] "Energy distance", https://en.wikipedia.org/wiki/Energy_distance
    [2] Szekely "E-statistics: The energy of statistical samples." Bowling
        Green State University, Department of Mathematics and Statistics,
        Technical Report 02-16 (2002).
    [3] Rizzo, Szekely "Energy distance." Wiley Interdisciplinary Reviews:
        Computational Statistics, 8(1):27-38 (2015).
    [4] Bellemare, Danihelka, Dabney, Mohamed, Lakshminarayanan, Hoyer,
        Munos "The Cramer Distance as a Solution to Biased Wasserstein
        Gradients" (2017). :arXiv:`1705.10743`.

"""
function cdf_distance(u::ECDF, v::ECDF, p::Int=1)

# values are pre-sorted ascending
    uv = u.sorted_values
    vv = v.sorted_values
    # merge and sort all grid values
    # We insert a zero to ensure that the first element's
    # grid distance survives the diff operation
    uplusv = vcat(0., uv, vv)
    sort!(uplusv)
    @debug "uplusv" uplusv
    Δx = diff(uplusv, dims = 1)
    @debug "Delta Force!" Δx

# Shortcut return in case the grids are identical
    if uv == vv # we can safely remove zeros from the delta
        Δx = filter(v->!iszero(v), Δx)
        UA = u(uv) # calculate u,v CDFs on their separate grids
        VA = v(vv)
    else
        A = uplusv[1:end-1]
        UA = u(A) # Calculate CDFs on whole domain
        VA = v(A)
    end

    if p == 1
        return sum(abs.(UA .- VA) .* Δx)
    elseif p == 2
        return sqrt(sum( (UA .- VA) .^ 2 .* Δx ))
    else
        return sum(abs.((UA .- VA)) .^ p .* Δx) ^ 1/p
    end
end

# *********** CONSTANTS ************

"""
Constant to use when specifying the Square Euclidean or L2 distance metric for a sequencing run.

    The L2 metric measures the sum of squares of the Euclidean or Manhattan-style 
    distance between two vectors of values. Because it operates only on values and not 
    on their underlying grid or axis, the L2 metric is less sensitive to the positions of  
    values along their grid than the EMD and Energy metrics.
    ``
        L2(x,y) = Σ ||x - y||²
    ``
See [`Distances.SqEuclidean`](@Ref)

"""
const L2 = Distances.SqEuclidean()

"""
Kullbach-Leibler Divergence metric.
"""
const KLD = Distances.KLDivergence()

"""
Monge-Wasserstein a.k.a. 1-p Wasserstein a.k.a. Earth Mover's Distance (EMD) metric. Sensitive to underlying 
grid. Default is unit grid taken from the axes of the data vector.
"""
const WASS1D = EMD()

"""
#TODO Fix this comment. Energy metric as defined by Szekely. (CHECK)

[wikipedia article](https://en.wikipedia.org/wiki/Energy_distance)
"""
const ENERGY = Energy()

"Convenience constant to represent L2, WASS1D, KLD, and ENERGY."
const ALL_METRICS = (L2, WASS1D, KLD, ENERGY)
