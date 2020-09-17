# SequencerJ.jl

```@meta
CurrentModule = SequencerJ
```

## Installation

SequencerJ is *not yet* a registered package, but installing it requires mere whispers of additional keystroke beyond those for a regular package:

```julia

`]` add "https://github.com/turingtest37/SequencerJ.jl/"

```
The `]` key switches the REPL into package manager mode. The `add` command accepts a URL. The  `delete` key gets you back to the REPL.

Alternatively you can stay in the REPL:

```@repl
using Pkg; Pkg.add(PackageSpec(url="https://github.com/turingtest37/SequencerJ.jl/"));
```
You may get **WARN**INGs upon compilation. You can safely ignore them for most purposes, but if you are developing SequencerJ locally and use the `Revise` package, note that you may have to restart your Julia environment more often than usual.

## Using SequencerJ

Getting started with SequencerJ is straightforward. First, we need to wrangle some data to analyze. For a quick (but unpromising) start, let's use a random array.
```@repl 1
using SequencerJ

A = rand(100,50);
r = sequence(A);
```
The `sequence` method is the primary interface to SequencerJ. Only the input matrix or vector A is required. All other options use their defaults: `scales=nothing` enables autoscaling to find the best row segmentation; `metrics=ALL_METRICS` enables all supported measurement; `grid=nothing` forces automatic grid creation.
See [`sequence`](@ref)

Data and statistics from a sequencing run are encapsulated in a [`SequencerResult`](@ref).

```@repl 1
r
```

Use SequencerJ's accessor functions to get the details of a run. The result data sequence (1-based *column* indices, not 0-based row indices as in Sequencer for python) can be obtained with the [`order`](@ref) function.

```@repl 1

order(r)

```

A key feature of the Sequencer algorithm is the calculated *elongation* of the minimum spanning tree that describes the best fit sequence. See [`elong`](@ref)

```@repl 1
elong(r)
```
### A Better Example
Random data looks pretty boring from every direction, so let's apply the Sequencer to something more enjoyable: trying to put back together an image whose columns (of pixels) have been shuffled. This trick will even impress your mom!

The `resources` folder contains test images you can use to play with SequencerJ.

```@example 2
using SequencerJ
using Random
using Images

img = load(joinpath(@__DIR__, "..", "..","resources","bread.jpeg"))
```

Color is nice but grayscale is easier.

```@example 2
imgg = (Gray.(img))
```

Let's slice up this loaf of bread with a different kind of knife...

```@example 2
A = convert(Matrix{Float32}, imgg);
shuff_idx = shuffle(axes(A,2));
prettyp(shuff_idx, 3)
```

Reordering the image with a shuffled column index creates a lot of bread crumbs!

```@example 2
As = A[:, shuff_idx]
colorview(Gray, As)
```


Now let's put back together the pieces. We run SequencerJ against the shuffled matrix,
and gather the best sequence as a column index.

```@example 2
seqres = sequence(As, metrics=(KLD, L2))
bestguessidx = SequencerJ.order(seqres)
```

Drum roll, please!

```@example 2
colorview(Gray, As[:, bestguessidx])
```

Oops, the image is mirrored. That's easy to fix:

```@example 2
colorview(Gray, As[:, reverse(bestguessidx)])
```


Voilà!

!!! warning "Failure is an option"
    The Sequencer sometimes fails to reassemble images completely. Usually this is because an image is too small along the row dimension and therefore does not contain enough data to allow the algorithm to discriminate correctly between similar columns of data. To ensure best performance, images should be presented in "portrait" mode, with more rows than columns.
    
See a full example of portrait mode here [Data in Portrait Mode](@ref)


## Metrics

To make sense of A, we must choose which statistical distance metrics to apply. Each metric compares A pairwise, column by column, to create a distance matrix. The matrix is analyzed using graph theory to identify an *optimal* ordering of columns in A.

This optimal sequence represents the least-cost path through the distance matrix and hence the closest affinities between columns, as proximity is the inverse of distance. See [The Sequencer Algorithm](@ref)

SequencerJ currently supports four metrics:

  * L2 - [`Euclidian`](@ref)
  * Earth Mover's Distance - [`EMD`](@ref)
  * Kullback-Lubler divergence - [`KLDivergence`](@ref)
  * Szekely's energy distance - [`Energy`](@ref)


With no other arguments, `sequence` applies all algorithms it knows to the data (`metrics = ALL_METRICS`). Distance algorithms may be specified individually or in groups with the `metrics` keyword.

```julia
r = sequence(A, scales=(1,3), metrics=(L2,ENERGY))
```

To specify a single metric, write it as a 1-tuple:

```julia
r = sequence(A, scales=(1,3), metrics=(L2,))
```

## Scales

The default `scale` is chosen by an "autoscaling" algorithm that samples the columns, runs the Sequencer on this subset at several scales (currently, the Fibonacci sequence), and compares them. The scale producing the greatest elongation is chosen to run on the full data set. You may force autoscaling with (`scales=nothing`).

```julia
sequence(A, scales=nothing)
```

Scale means the number of parts into which the data is partitioned. Each section or *chunk* contains approximately `size(A,1)/scale` elements. For example, 100 rows at scale 3 will result in chunks of 33, 33, and 34 rows. Set scales manually using the `scales` keyword:

```julia
r = sequence(A, scales=(1,3))
```

To choose only one scale, write it as a 1-tuple:

```julia
r = sequence(A, scales=(4,))
```

## Output

Use accessor functions to get details from the `SequencerResult`.

```@repl accessors
using SequencerJ #hide
A = rand(2000, 100);
r = sequence(A);
```

Best column sequence: [`order`](@ref)
```@repl accessors
order(r)
```

Final elongation: [`elong`](@ref)
```@repl accessors
elong(r)
```

Final elongation: [`mst`](@ref)
```@repl accessors
mst(r)
```

Final distance matrix: [`D`](@ref)
```@repl accessors
D(r)
```

Per-segment intermediate results:
```@repl accessors
using Base.Iterators
segDict = r.EOSeg;
for (k,l) in Iterators.take(keys(segDict), 3)
    η, mst = segDict[(k,l)]
end
segDict[KLD]
η, mst = segDict[(KLD,2)]
η = first(segDict[(KLD,2)])
```

Collect the mean elongations across segments for each metric+scale
```@repl accessors
collect(StatsBase.mean(first(v)) for (k,v) in r)
```

In a similar fashion, get final elongations and the MST for each metric+scale:
```@repl accessors
rk = r.EOAlgScale
```

## The Sequencer Algorithm

The algorithm is described in detail [in the paper available on Arxiv](https://arxiv.org/abs/2006.13948)

    @misc{baron2020extracting,
    title={Extracting the main trend in a dataset: the Sequencer algorithm},
    author={Dalya Baron and Brice Ménard},
    year={2020},
    eprint={2006.13948},
    archivePrefix={arXiv},
    primaryClass={cs.LG},
    year=2020
}

As another example of how the Sequencer algorithm may be used, see [Mapping Earth's deepest secrets](https://science.sciencemag.org/content/368/6496/1183).
