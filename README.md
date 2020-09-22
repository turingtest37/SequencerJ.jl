# SequencerJ

__SequencerJ__ is a pure [Julia](https://julialang.org/) implementation of the [Sequencer Algorithm](https://github.com/dalya/Sequencer/), a data analysis tool to identify and extract the principal trends in a set of 1-d data vectors.

Getting started with SequencerJ is easy. From the Julia REPL:
```julia
    julia> using Pkg; Pkg.add(PackageSpec(url="https://github.com/turingtest37/SequencerJ.jl/"))
    [...]
    julia> using SequencerJ
    [ Info: Precompiling SequencerJ [348581b9-6e84-42e0-ac4e-fe9177c221e6]
    [...]
```
You may get **WARN**INGs upon compilation. You can safely ignore them for most purposes, but if you are developing SequencerJ locally and use the `Revise` package, note that you may have to restart your Julia environment more often than usual.

```julia
    julia> A = rand(50,100); #some data to process. 

    julia> m = ALL_METRICS
    (Euclidean(0.0), EMD(nothing), KLDivergence(), Energy(nothing))

    julia> s = (1,2,4)
    (1, 2, 4)

    julia> seqres = sequence(A; metrics=m, scales=s)
    ┌ Info: Sequencing data with
    │     shape: (50, 100)
    │     metric(s): (Euclidean(0.0), EMD(nothing), KLDivergence(), Energy(nothing))
    └     scale(s): (1, 2, 4)
    [...]
```

`seqres` is a `SequencerResult` type that may be used to retrieve results of the sequencing run. The `order` function returns the best reordering of the data columns that was found.

```julia
    julia> bestseq = order(seqres)
    100-element Array{Int64,1}:
    10
    15
    13
    [...]
```

The Sequencer also calculates a fitness coefficient `η` that can be used to compare quality of solutions using various metrics and scales against the same data. Bigger is better. η is returned by the `elong` function.
```julia
    julia> eta = elong(seqres)
    6.2345
```


The paper that describes the Sequencer algorithm and its applications can be found 
on Arxiv: [https://arxiv.org/abs/2006.13948].
