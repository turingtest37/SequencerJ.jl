# SequencerJulia

SequencerJulia is a pure Julia implementation of the Sequencer Algorithm, an analysis tool to identify and extract the principal trends in a set of 1-d data vectors.

Using SequencerJulia is easy:

```
jldoctest

julia> using Pkg; Pkg.add(PackageSpec(url="https://github.com/turingtest37/SequencerJulia.jl/"))
[...]

julia> using SequencerJulia
[ Info: Precompiling SequencerJulia [348581b9-6e84-42e0-ac4e-fe9177c221e6]
[...]

julia> A = rand(50, 100)
50×100 Array{Float64,2}:
[...]

julia> m = ALL_METRICS
(SqEuclidean(0.0), EMD(nothing), KLDivergence(), Energy(nothing))

julia> s = (1,2,4)
(1, 2, 4)

julia> sequence(A; metrics=m, scales=s)
┌ Info: Sequencing data with
│     shape: (50, 100)
│     metric(s): (SqEuclidean(0.0), EMD(nothing), KLDivergence(), Energy(nothing))
└     scale(s): (1, 2, 4)
[ Info: SqEuclidean(0.0) at scale 1: η = 1.338
[ Info: SqEuclidean(0.0) at scale 2: η = 1.702
[ Info: SqEuclidean(0.0) at scale 4: η = 2.521
[ Info: EMD(nothing) at scale 1: η = 28.43
[ Info: EMD(nothing) at scale 2: η = 24.65
[ Info: EMD(nothing) at scale 4: η = 22.22
[ Info: KLDivergence() at scale 1: η = 19.65
[ Info: KLDivergence() at scale 2: η = 15.43
[ Info: KLDivergence() at scale 4: η = 15.98
[ Info: Energy(nothing) at scale 1: η = 7.413
[ Info: Energy(nothing) at scale 2: η = 8.564
[ Info: Energy(nothing) at scale 4: η = 10.57
[ Info: Final average elongation: 1.361e+37
Sequencer Result: η = 1.361e+37, order = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100] 

```


The paper that describes the Sequencer algorithm and its applications can be found 
on Arxiv: [https://arxiv.org/abs/2006.13948].
```
@misc{baron2020extracting,
    title={Extracting the main trend in a dataset: the Sequencer algorithm},
    author={Dalya Baron and Brice Ménard},
    year={2020},
    eprint={2006.13948},
    archivePrefix={arXiv},
    primaryClass={cs.LG},
    year=2020
}
```