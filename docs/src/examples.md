# Examples of SequencerJ




## Data in Portrait Mode

SequencerJ works best with data going down the columns.

```@example

using SequencerJ
using Random
using Images

img = Gray.(load(joinpath(@__DIR__, "..","..","resources","bread.jpeg")))
A = convert(Matrix{Float32}, img)
M,N = size(A)
A = M < N ? permutedims(A) : A
shuff_idx = shuffle(axes(A,2))
As = A[:, shuff_idx]
seqres = sequence(As, metrics=(KLD, L2)); #autoscale
colorview(Gray, permutedims(As[:, order(seqres)])) # re-permute the matrix to get back landscape

```