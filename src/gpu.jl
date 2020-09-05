const USECUDA = CuArrays.functional()
ENV["JULIA_CUDA_VERBOSE"] = USECUDA
ENV["JULIA_CUDA_SILENT"] = !USECUDA

allowscalar(false)

wrap(A::AbstractArray) = (USECUDA && !(typeof(A) <: CuArray)) ? CuArray(A) : A
