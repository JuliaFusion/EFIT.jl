# EFIT

[EFIT (Equilibrium Fitting)](https://fusion.gat.com/theory/Efit) is a computer code developed to translate measurements from plasma diagnostics into useful information like plasma geometry, stored energy, and current profiles.


EFIT.jl provides basic functionality for reading [EFIT GEQDSK](https://fusion.gat.com/theory/Efitgeqdsk) files.

```
julia> using EFIT

julia> g = readg(EFIT.test_gfile);

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

julia> triangularity(g)
(0.059614676027860296, 0.05822145848512557)

julia> ellipticity(g)
1.475784591289634

julia> elongation(g)
1.475784591289634

julia> major_radius(g)
1.648852555

julia> minor_radius(g)
0.633397135

julia> aspect_ratio(g)
2.6031891587258285
```
