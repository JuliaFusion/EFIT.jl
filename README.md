[![Build Status](https://travis-ci.org/JuliaFusion/EFIT.jl.svg?branch=master)](https://travis-ci.org/JuliaFusion/EFIT.jl)

# EFIT

[EFIT (Equilibrium Fitting)](https://fusion.gat.com/theory/Efit) is a computer code developed to translate measurements from plasma diagnostics into useful information like plasma geometry, stored energy, and current profiles.


`EFIT.jl` provides basic functionality for reading [EFIT GEQDSK](https://fusion.gat.com/theory/Efitgeqdsk) files.

```julia-repl
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


`EFIT.jl` can convert `dd` into GEQDSK files as in the following example:
```julia
import EFIT
import EFIT.IMASdd

# Load example dd
filename = joinpath(dirname(dirname(pathof(EFIT))), "test", "test_dd_eq.json")
dd = IMASdd.json2imas(filename);

# First, convert equilibrium time_slices in dd to a list of `GEQDSK` struct
ggs = EFIT.imas2geqdsk(dd)

# Then, write into files
for gg in ggs
    # default file name has "g00000.xxxxx" format, where xxxxx is time in `ms`
    file_name = gg.file
    EFIT.writeg(gg, file_name)
end
```