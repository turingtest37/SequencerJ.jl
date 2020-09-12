# Metrics

```@meta
CurrentModule = SequencerJ
```

## Types
SequencerJ currently supports the four algorithms originally cited in the article and supported by the `python` version of the Sequencer. These are:



```@docs
SqEuclidean
KLDivergence
EMD
Energy
```

Each algorithm has associated to it a constant that represents a usable default type instance (with a default grid):

## Constants

```@docs
L2
KLD
WASS1D
ENERGY
ALL_METRICS
```

## Usage

### Data Orientation
SequencerJ operates on pairwise *columns* of data in the given matrix (or vector of vectors). It  builds an N x N distance matrix for further analysis, given an N-column input. It is therefore important to properly orient the data matrix so that M observations run downward along N columns, with each column representing a separate data set or series.

This can be performed efficiently with `permutedims`:
```@example
A = rand(10,5)
M,N = size(A)
A = M < N ? permutedims(A) : A
```

While a 2-D data set such as a photograph may be shuffled and reassembled in either orientation,
the Sequencer algorithm will work better with more rows and fewer columns (i.e. M > N).

### Autoscaling

The Sequencer algorithm's processing time rises as N x N, posing a challenge for data exploration with larger data sets. It may require considerable experimentation with multiple scales and metrics before finding an optimum combination. To assist with finding an appropriate scale, SequencerJ offers an "autoscaling" option. This is, in fact, the default behavior with the keyword argument `scales=nothing`.

In autoscaling, SequencerJ first performs the sequencing function at several scales on a subset of columns. The scale producing the greatest elongation (see [`elong`](@ref)) is then used on the full dataset. This is a feature specific to the Julia port of the Sequencer algorithm.

At present, there is no "autometric" option to choose the best metric.

### Metrics and Grids

Setting the metrics keyword option as
```julia
sequence(A, metrics=(L2, ENERGY))
```
means to use `Distances.SqEuclidean()` and `SequencerJ.Energy()`.


With equal effect, you can write:
```julia
sequence(A, metrics=(SqEuclidean(), Energy())
```


You can specify a threshold for the `SqEuclidean` calculations. To use, e.g., `1e-7` means round to zero any distance calculation results smaller than `1e-7`.
```julia
sequence(A, metrics=(SqEuclidean(1e-7), Energy())
```


To specify a non-unit data grid, set `grid` as a `Vector{AbstractFloat}` in the range [0.0-1.0) and use the Distance type, not the constant.
```example grid
using SequencerJ #hide
A = rand(100,100);
M,N = size(A)
g = collect(0.0:0.5:N-1)
sequence(A, metrics=(Energy()), grid=g)
```

To these four metrics, the intention is to add spectral methods for complex (i.e. 2-D) inputs and any other Distances that make sense for sequencing 1-D or 2-D data. I would love community help with this!
