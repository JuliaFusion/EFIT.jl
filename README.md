# EFIT

[EFIT (Equilibrium Fitting)](https://fusion.gat.com/theory/Efit) is a computer code developed to translate measurements from plasma diagnostics into useful information like plasma geometry, stored energy, and current profiles.


EFIT.jl provides basic functionality for reading [EFIT GEQDSK](https://fusion.gat.com/theory/Efitgeqdsk) files.

```
julia> using EFIT

julia> g = readg("test/g000001.01000")
GEQDSKFile: "test/g000001.01000"

julia> g.fpol
101-element Array{Float64,1}:
 -3.38249
 -3.38691
 -3.39044
 -3.3932
 -3.39546
 -3.39738
 -3.39903
 ...
```
