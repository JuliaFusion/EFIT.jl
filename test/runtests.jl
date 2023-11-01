using Test
using EFIT

g = readg(EFIT.test_gfile)
g2 = readg(EFIT.test_gfile2)

lt, ut = triangularity(g)
@test lt == 0.059614676027860296
@test ut == 0.05822145848512557

@test ellipticity(g) == 1.475784591289634
@test major_radius(g) == 1.648852555
@test minor_radius(g) == 0.633397135
@test aspect_ratio(g) == 2.6031891587258285

@testset "find xpoint" begin
    psin_tolerance = 0.002
    xrs, xzs, xpsins, xseps = EFIT.x_points(g, psin_tolerance=psin_tolerance)
    @test length(xrs) == 0  # this file is limited, so no x-points should be found
    @test length(xrs) == length(xzs) == length(xpsins) == length(xseps)

    xrs, xzs, xpsins, xseps = EFIT.x_points(g2, psin_tolerance=psin_tolerance)
    @test length(xrs) > 0
    @test length(xrs) == length(xzs) == length(xpsins) == length(xseps)
    @test xseps[1] == 1
    @test abs(xpsins[1] - 1.0) < psin_tolerance
    @test xrs[1] > 0
    @test (xrs[1] > 1.25) & (xrs[1] < 1.26)
    @test (xzs[1] > -1.17) & (xzs[1] < -1.15)

    @test xseps[2] == 2
    @test xpsins[2] > (1 + psin_tolerance)
    @test (xrs[2] > 1.2) & (xrs[2] < 1.33)
    @test (xzs[2] > 1.04) & (xzs[2] < 1.17)
end
