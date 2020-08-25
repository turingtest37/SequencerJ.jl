# Metrics

```@meta
CurrentModule = SequencerJulia
```

## EMD
```@docs
EMD()
EMD(u,v)
EMD(u,v,uw,vw)
(m::EMD)(u,v)
(m::EMD)(u,v,uw,vw)
emd(u,v)
emd(u,v,uw,vw)
```

## Energy
```@docs
Energy()
Energy(u,v)
Energy(u,v,uw,vw)
(m::Energy)(u,v)
(m::Energy)(u,v,uw,vw)
energy(u::AbstractVector, v=u)
energy(u,v,uw,vw)

```

### L2
```@docs
Distances.L2
```

### KL-Divergence
```@docs
Distances.KLDivergence
```
