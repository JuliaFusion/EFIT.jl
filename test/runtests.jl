using Test
using EFIT

g = readg(EFIT.test_gfile)

lt, ut = triangularity(g)
@test lt == 0.059614676027860296
@test ut == 0.05822145848512557

@test ellipticity(g) == 1.475784591289634
@test major_radius(g) == 1.648852555
@test minor_radius(g) == 0.633397135
@test aspect_ratio(g) == 2.6031891587258285

