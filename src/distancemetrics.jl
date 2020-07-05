
struct EMD{T <: Union{Nothing,AbstractVector}} <: SemiMetric
    grid::T
end

function EMD(u,v)
    e = EMD(Array(first(axes(u))))
    length(u) ==  length(v) || error("u and v must have same length. Got u $(length(u)) and v $(length(v))")
    # size(u) == size(v) || error("Only 1-d data vectors are supported.")
    wasserstein1d(ecdf(e.grid; weights = u), ecdf(e.grid; weights = v))
end

function EMD(u,v,uw,vw)
    length(u) ==  length(uw) || error("u and u weights must have same length. Got u $(length(u)) and uw $(length(uw))")
    length(v) ==  length(vw) || error("v and v weights must have same length. Got v $(length(v)) and vw $(length(vw))")
    # size(u) == size(v) || error("Only 1-d data vectors are supported.")
    wasserstein1d(ecdf(u; weights = uw), ecdf(v; weights = vw))
end

emd(a::AbstractArray, b::AbstractArray; grid=nothing) = EMD(grid)(a, b)

"""
Returns the Euclidean distance between all pairs of
elements in the given list.

DM = Σ ((X[i] - X[j])^2) for i,j=i+1 indexes on X

🚟
"""
function l2(objs; parallel = false)
    A = hcat(objs...)'
    return pairwise(SqEuclidean(), A, dims=1)
end

"""
Calculates the KL Divergence between pairs of elements
in the given list and returns them as a symmetric distance matrix.

This function assumes as an argument a vector of arrays.

"""
function kl(objs)
    A = hcat(objs...)'
    @debug "A" A
    return pairwise(KLDivergence(), A, dims=1)
end

function kl(u, v)
    A = hcat(u, v)
    @debug "A" A
    return pairwise(KLDivergence(), A, dims=1)
end

kl(A::AbstractMatrix) = pairwise(KLDivergence(), A, dims=1)

wasserstein1d(u, v, uw, vw) = wasserstein1d(ecdf(u; weights=uw), ecdf(v; weights=vw))

wasserstein1d(u::ECDF, v::ECDF) = cdf_distance(u, v, 1)

energy(u::ECDF, v::ECDF) = sqrt(2) * cdf_distance(u, v, 2)

struct Energy{T <: Union{Nothing,AbstractVector}} <: SemiMetric
    grid::T
end

function Energy(u,v)
    e = Energy(Array(first(axes(u))))
    length(u) ==  length(v) || error("u and v must have same length. Got u $(length(u)) and v $(length(v))")
    # size(u) == size(v) || error("Only 1-d data vectors are supported.")
    energy(ecdf(e.grid; weights = u), ecdf(e.grid; weights = v))
end

"""
Returns the statistical distance between two cumulative
distributions u et v. With dimension parameter p = 1,
the distance is the Wasserstein Distance. With p=2,

"""
function cdf_distance(u::ECDF, v::ECDF, p::Int=1)

# values are pre-sorted ascending
    uv = Float64.(u.sorted_values)
    vv = Float64.(v.sorted_values)
# weights are pre-sorted in ECDF to match values
    uw = Float64.(u.weights)
    vw = Float64.(v.weights)

# merge and sort all values
    uplusv = vcat(uv, vv)
    sort!(uplusv)

# the underlying grid space on which we operate
    deltas = diff(uplusv, dims = 1)
# creates a n-1 size result for a length n list
    @debug "Delta Force!" deltas

# cumulated sum of weights => CDF
    sumuw = vcat([0.], cumsum(uw))
    sumvw = vcat([0.], cumsum(vw))

# add 1 to each index to accomodate [0] added to cumulative sums above
    ucdf_i = collect(searchsortedlast(uv, x) + 1 for x in uplusv[1:end-1])
    vcdf_i = collect(searchsortedlast(vv, x) + 1 for x in uplusv[1:end-1])
# normalize CDF to range [0,1]
    ucdf = sumuw[ucdf_i] ./ last(sumuw)
    vcdf = sumvw[vcdf_i] ./ last(sumvw)
    @debug "ucdf" ucdf
    @debug "vcdf" vcdf

    if p == 1
        return sum(abs.(ucdf .- vcdf) .* deltas)
    elseif p == 2
        return sqrt(sum((ucdf .- vcdf) .^ 2 .* deltas))
    else
        return sum(abs.(ucdf .- vcdf) .^ p .* deltas) ^ 1/p
    end
end


"""
Calculates the 1-d Earth Mover Distance between
elements of the given list. Uses the 1-d Wasserstein Distance.
"""
function emd(weights::AbstractArray, grid::Union{Nothing,AbstractArray} = nothing)
    DM = zeros(Float64,length(weights),length(weights))

    grid = isnothing(grid) ? collect(i for i in first(axes(first(weights)))) : grid
    @debug "Grid" grid

    size(grid) == size(first(weights)) || throw(ArgumentError("weights is messed up."))

    collect(DM[i,j] = DM[j,i] = wasserstein1d(ecdf(grid; weights=weights[i]), ecdf(grid; weights=weights[j])) for i in 1:length(weights) for j in (i+1):length(weights));
    return DM

end
