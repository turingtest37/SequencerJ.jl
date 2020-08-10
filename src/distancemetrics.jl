
struct EMD <: SemiMetric end

"""
Returns fhe Earth Mover Distance (EMD) a.k.a the 1-Wasserstein distance
between the two given vectors (which are treated as weights on a grid
defined by the first axis of u).
"""
function EMD(u, v = u)
    grid = Array(first(axes(u)))
    EMD(grid, grid, u, v)
end

function EMD(u,v,uw,vw)
    ndims(u) == 1 && ndims(u) == 1 || error("Only 1-d data vectors are supported.")
    length(u) == length(uw) || error("u and u weights must have same length. Got u $(length(u)) and uw $(length(uw))")
    length(v) == length(vw) || error("v and v weights must have same length. Got v $(length(v)) and vw $(length(vw))")
    cdf_distance(ecdf(u; weights = uw), ecdf(v; weights = vw), 1)
end

(::EMD)(u,v=u) = EMD(u,v)
(::EMD)(u,v,uw,vw) = EMD(u,v,uw,vw)

emd(u::AbstractArray, v = u; grid=nothing) = EMD(u, v)

emd(u,v,uw,vw) = EMD(u,v,uw,vw)


# ******** ENERGY *********
# To get distance matrix, use pairwise(Energy(), A, dims=[1|2])
# A has to be a Matrix. Use dims to indicate direction of columns
# A = hcat(objs...)'
# @debug "A" A
# return pairwise(Energy(), A, dims=1)

struct Energy <: SemiMetric end

# enables dispatching on the type directly
(::Energy)(u,v) = Energy(u,v)
(::Energy)(u,v,uw,vw) = Energy(u,v,uw,vw)

function Energy(u,v,uw,vw)
    length(u) ==  length(uw) || error("u and u weights must have same length. Got u $(length(u)) and uw $(length(uw))")
    length(v) ==  length(vw) || error("v and v weights must have same length. Got v $(length(v)) and vw $(length(vw))")
    ndims(u) == 1 && ndims(v) == 1 || error("Only 1-d data vectors are supported.")
    sqrt(2) * cdf_distance(ecdf(u; weights = uw), ecdf(v; weights = vw), 2)
end

function Energy(u,v)
    gridu = Array(first(axes(u)))
    gridv = Array(first(axes(v)))
    Energy(gridu,gridv,u,v)
end

energy(u, v, uw, vw) = Energy(u, v, uw, vw)

energy(u::AbstractVector, v=u) = Energy(u,v)

# TODO Upgrade to EmpiricalCDF package instead of StatsBase.ECDF
"""
Returns the statistical distance between two cumulative
distributions u et v. With dimension parameter p = 1,
the distance is the 1 Wasserstein (or Earth Mover's) Distance.
With p=2, the distance is Energy.

"""
function cdf_distance(u::ECDF, v::ECDF, p::Int=1)

# values are pre-sorted ascending
    uv = Float32.(u.sorted_values)
    vv = Float32.(v.sorted_values)
# weights are pre-sorted in ECDF to match values
    uw = Float32.(u.weights)
    vw = Float32.(v.weights)

# merge and sort all values
    uplusv = vcat(uv, vv)
    sort!(uplusv)
    @debug "uplusv" uplusv
# the underlying grid space on which we operate
    deltas = diff(uplusv, dims = 1)
# creates a n-1 size result for a length n list
    @debug "Delta Force!" deltas

    u2 = ecdf(push!(uv, 0.0); weights=push!(uw, 0.0))
    v2 = ecdf(push!(vv, 0.0); weights=push!(vw, 0.0))
    A = uplusv[1:end-1]
    @debug "A = uplusv[1:end-1]" A

    if p == 1
        # return sum(abs.(ucdf .- vcdf) .* deltas)
        return sum(abs.(u2(A) .- v2(A)) .* deltas)
    elseif p == 2
        return sqrt(sum( (u2(A) .- v2(A)) .^ 2 .* deltas ))
        # return sqrt(sum((ucdf .- vcdf) .^ 2 .* deltas))
    else
        return sum(abs.(u.(A) - v.(A)) .^ p .* deltas) ^ 1/p
    end
end
