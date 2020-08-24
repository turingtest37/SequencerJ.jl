# SequencerJulia

SequencerJulia is a pure [Julia](https://julialang.org/) implementation of the [Sequencer Algorithm](https://github.com/dalya/Sequencer/), a data analysis tool to identify and extract the principal trends in a set of 1-d data vectors.

Using SequencerJulia is straightforward. From the Julia REPL:

    julia> using Pkg; Pkg.add(PackageSpec(url="https://github.com/turingtest37/SequencerJulia.jl/"))
    [...]

    julia> using SequencerJulia
    [ Info: Precompiling SequencerJulia [348581b9-6e84-42e0-ac4e-fe9177c221e6]
    [...]
    # You may get WARNINGs upon compilation. You can safely ignore them.

    julia> A = rand(50, 100)
    50×100 Array{Float64,2}:
    [...]

 Measure A using L2, Earth Mover's Distance, the Kullback-Lubler divergence, and Szekely's energy metric.
 See [`Distances.SqEuclidian`](@ref), [`EMD`](@ref), [`Distances.KLDivergence`](@ref), [`Energy`](@ref) 
    julia> m = ALL_METRICS
    (SqEuclidean(0.0), EMD(nothing), KLDivergence(), Energy(nothing))

    julia> s = (1,2,4)
    (1, 2, 4)

    julia> sequence(A; metrics=m, scales=s)
    ┌ Info: Sequencing data with
    │     shape: (50, 100)
    │     metric(s): (SqEuclidean(0.0), EMD(nothing), KLDivergence(), Energy(nothing))
    └     scale(s): (1, 2, 4)
    [ Info: SqEuclidean(0.0) at scale 1: η = 2.716 (3.7s)
    [ Info: SqEuclidean(0.0) at scale 2: η = 3.867 (0.012s)
    [ Info: SqEuclidean(0.0) at scale 4: η = 3.528 (0.038s)
    [ Info: EMD(nothing) at scale 1: η = 6.83 (1.4s)
    [ Info: EMD(nothing) at scale 2: η = 5.385 (0.12s)
    [ Info: EMD(nothing) at scale 4: η = 4.286 (0.19s)
    [ Info: KLDivergence() at scale 1: η = 3.355 (0.17s)
    [ Info: KLDivergence() at scale 2: η = 3.745 (0.021s)
    [ Info: KLDivergence() at scale 4: η = 2.867 (0.052s)
    [ Info: Energy(nothing) at scale 1: η = 6.086 (0.19s)
    [ Info: Energy(nothing) at scale 2: η = 5.036 (0.11s)
    [ Info: Energy(nothing) at scale 4: η = 4.665 (0.2s)
    [ Info: Final average elongation: 3.036
    [ Info: Final ordering: [76, 47, 87, 56, 14, 7, 68, 93, 88, 73, 71, 66, 95, 72, 34, 30, 53, 11, 4, 99, 3, 92, 78, 69, 67, 90, 39, 98, 80, 77, 42, 62, 59, 48, 45, 31, 15, 89, 22, 51, 94, 64, 41, 28, 58, 91, 82, 52, 27, 44, 10, 35, 50, 24, 2, 75, 96, 63, 55, 33, 26, 43, 1, 12, 18, 17, 5, 60, 49, 54, 85, 84, 81, 25, 40, 32, 6, 46, 100, 38, 19, 37, 36, 20, 16, 23, 57, 13, 97, 61, 9, 29, 8, 86, 74, 70, 65, 21, 79, 83]
    Sequencer Result: η = 3.036, order = [76, 47, 87, 56, 14, 7, 68, 93, 88, 73, 71, 66, 95, 72, 34, 30, 53, 11, 4, 99, 3, 92, 78, 69, 67, 90, 39, 98, 80, 77, 42, 62, 59, 48, 45, 31, 15, 89, 22, 51, 94, 64, 41, 28, 58, 91, 82, 52, 27, 44, 10, 35, 50, 24, 2, 75, 96, 63, 55, 33, 26, 43, 1, 12, 18, 17, 5, 60, 49, 54, 85, 84, 81, 25, 40, 32, 6, 46, 100, 38, 19, 37, 36, 20, 16, 23, 57, 13, 97, 61, 9, 29, 8, 86, 74, 70, 65, 21, 79, 83]


# Use accessor functions to get details of the results

Elongation coefficent of the final MST

    julia> elong(result)
    3.0355999999999996

Final vertex sequence (*column* indices, not row as in Sequencer for python)

    julia> order(result)
    100-element Array{Int64,1}:
    76
    47
    87
    [...]
    79
    83

 Final distance matrix   

    julia> D(result)
    100×100 SparseArrays.SparseMatrixCSC{Float64,Int64} with 10000 stored entries:
    [1  ,   1]  =  1.0e-6
    [2  ,   1]  =  Inf
    [...]
    [99 , 100]  =  Inf
    [100, 100]  =  1.0e-6

Per-segment intermediate results:

    julia> r = result.EOSeg;

    julia> r[KLD,2]
    (Any[4.139200000000001, 2.7687999999999997], Any[{100, 99} directed simple Int64 graph, {100, 99} directed simple Int64 graph])

    julia> η, mst = r[KLD,2]
    (Any[4.139200000000001, 2.7687999999999997], Any[{100, 99} directed simple Int64 graph, {100, 99} directed simple Int64 graph])

    julia> η = first(r[KLD,2])
    2-element Array{Any,1}:
    4.139200000000001
    2.7687999999999997

Collect the mean elongations across segments for each metric+scale

    julia> collect(StatsBase.mean(first(v)) for (k,v) in r)
    12-element Array{Float64,1}:
    6.0466
    8.0551
    7.6366
    [...]

In a similar fashion, get final elongations and the MST for each metric+scale

    julia> rk = result.EOAlgScale
    Dict{Tuple,Any} with 12 entries:
    (SqEuclidean(0.0), 4) => (3.5284, {100, 99} directed simple Int64 graph)
    (EMD(nothing), 4)     => (4.2864, {100, 99} directed simple Int64 graph)
    (Energy(nothing), 4)  => (4.6652, {100, 99} directed simple Int64 graph)
    (SqEuclidean(0.0), 2) => (3.8672, {100, 99} directed simple Int64 graph)
    (KLDivergence(), 1)   => (3.3548, {100, 99} directed simple Int64 graph)
    (SqEuclidean(0.0), 1) => (2.716, {100, 99} directed simple Int64 graph)
    (EMD(nothing), 2)     => (5.3852, {100, 99} directed simple Int64 graph)
    (Energy(nothing), 1)  => (6.0862, {100, 99} directed simple Int64 graph)
    (KLDivergence(), 4)   => (2.8672, {100, 99} directed simple Int64 graph)
    (Energy(nothing), 2)  => (5.0356, {100, 99} directed simple Int64 graph)
    (EMD(nothing), 1)     => (6.83, {100, 99} directed simple Int64 graph)
    (KLDivergence(), 2)   => (3.745, {100, 99} directed simple Int64 graph)

More and better docs to follow!

The paper that describes the Sequencer algorithm and its applications can be found 
on Arxiv: [https://arxiv.org/abs/2006.13948].

@misc{baron2020extracting,
    title={Extracting the main trend in a dataset: the Sequencer algorithm},
    author={Dalya Baron and Brice Ménard},
    year={2020},
    eprint={2006.13948},
    archivePrefix={arXiv},
    primaryClass={cs.LG},
    year=2020
}
