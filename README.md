# SequencerJulia

SequencerJulia is a pure Julia implementation of the Sequencer Algorithm, an analysis tool to identify and extract the principal trends in a set of 1-d data vectors.

Using SequencerJulia is easy:

```
jldoctest

using Pkg; Pkg.add("https://github.com/turingtest37/SequencerJulia.jl/")

using SequencerJulia

A = rand(50, 100)
m = ALL_METRICS
s = (1,2,4)

sequencer(A; metrics=m, scales=s)

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