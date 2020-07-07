
struct EMD{T <: Union{Nothing,AbstractVector}} <: SemiMetric
    grid::T
end


"""
Returns fhe Earth Mover Distance (EMD) a.k.a the 1-Wasserstein distance
between the two given arrays (which are treated as weights on a grid
defined by the first axis of u).
"""
function EMD(u, v = u)
    grid = Array(first(axes(u)))
    EMD(grid, grid, u, v)
end

function EMD(u,v,uw,vw)
    length(u) ==  length(uw) || error("u and u weights must have same length. Got u $(length(u)) and uw $(length(uw))")
    length(v) ==  length(vw) || error("v and v weights must have same length. Got v $(length(v)) and vw $(length(vw))")
    length(size(u)) == 1 && length(size(u)) == 1 || error("Only 1-d data vectors are supported.")
    _wasserstein1d(ecdf(u; weights = uw), ecdf(v; weights = vw))
end

EMD() = EMD(nothing)

(::EMD)(u,v=u) = EMD(u,v)

emd(a::AbstractArray, b = a; grid=nothing) = EMD(grid)(a, b)

emd(u,v,uw,vw) = EMD(u,v,uw,vw)

"""
Returns a matrix of the squared euclidean distances between all pairs of
elements in the given list (typically a vector of vectors).

DM = Σ ((X[i] - X[j])^2) for i,j indexes on X

🚟
"""
function l2dm(objs)
    A = hcat(objs...)
    return pairwise(SqEuclidean(), A, dims=2)
end

l2(u,v) = SqEuclidean()(u,v)

"""
Calculates the KL Divergence between pairs of elements
in the given list and returns them as a symmetric distance matrix.

This function assumes as an argument a vector of vectors.

"""
function kl(objs::AbstractVector)
    A,B = objs
    @debug "A" A
    @debug "B" B
    return KLDivergence()(A,B)
end

function kldm(objs)
    A = hcat(objs...)
    @debug "A" A
    return pairwise(KLDivergence(), A, dims=2)
end

function kldm(u, v)
    A = hcat(u,v)
    @debug "A" A
    return pairwise(KLDivergence(), A, dims=2)
end

function kl(u, v)
    objs = hcat(u, v)
    @debug "objs" objs
    kl(objs)
    # return pairwise(KLDivergence(), A, dims=2)
end

kl(A::AbstractMatrix) = KLDivergence()(eachcol(A)...)

wasserstein1d(u, v, uw, vw) = _wasserstein1d(ecdf(u; weights=uw), ecdf(v; weights=vw))

_wasserstein1d(u::ECDF, v::ECDF) = cdf_distance(u, v, 1)


struct Energy{T <: Union{Nothing,AbstractVector}} <: SemiMetric
    grid::T
end

function Energy(u,v,uw,vw)
    length(u) ==  length(uw) || error("u and u weights must have same length. Got u $(length(u)) and uw $(length(uw))")
    length(v) ==  length(vw) || error("v and v weights must have same length. Got v $(length(v)) and vw $(length(vw))")
    length(size(u)) == 1 && length(size(u)) == 1 || error("Only 1-d data vectors are supported.")
    _energy(ecdf(u; weights = uw), ecdf(v; weights = vw))
end

function Energy(u,v)
    energy(u,v)
end

Energy() = Energy(nothing)

(::Energy)(u,v) = Energy(u,v)

_energy(u::ECDF, v::ECDF) = sqrt(2) * cdf_distance(u, v, 2)

function energy(u, v, uw, vw)
    _energy(ecdf(u; weights=uw), ecdf(v; weights=vw))
end

function energy(u::AbstractVector, v::AbstractVector)
    grid = Array(first(axes(u)))
    _energy(ecdf(grid; weights=u), ecdf(grid; weights=v))
end

function energy(objs)
    A = hcat(objs...)'
    @debug "A" A
    return pairwise(Energy(), A, dims=1)
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
