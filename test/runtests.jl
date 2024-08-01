using Test
using EFIT
import EFIT
import IMASdd

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
    # fluxinfo
    izmin = 5
    irmin = 4
    r = collect(LinRange(1.2, 1.4, 9))
    z = collect(LinRange(-1.45, -1.2, 9))
    rc = r[irmin]
    zc = z[izmin]
    lr = 0.0265
    lz = 0.0279
    psi = ((z .- zc) .^ 2 ./ (lz^2))' .- (r .- rc) .^ 2 ./ (lr^2) .+ 34.58
    br, bz, bpol, d = EFIT.fluxinfo(r, z, psi)
    @test size(br) == size(psi)
    minbp, minbploc = findmin(bpol)
    @test minbploc[1] == irmin
    @test minbploc[2] == izmin
    @test all(bpol .>= 0.0)

    # find_local_minimum
    nbh = 2
    bplm = EFIT.find_local_minimum(bpol, nbh)
    @test all(bplm[irmin-nbh:irmin+nbh, izmin-nbh:izmin+nbh] .== minbp)
    @test all(bplm[irmin+nbh+1:end, :] .!= minbp)
    @test all(bplm[1:irmin-nbh-1, :] .!= minbp)

    # x_point
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

@testset "parse header" begin
    headers = [
        "   EFITD    00/00/2008    #002296  0200ms Convert   0 257 513",
        "   EFITD   11/23/2020    #184833  3600             3  65  65",
    ]
    times = [
        0.2,
        3.6,
    ]
    expect_exception = [
        true,
        false,
    ]
    the_set_time = 0.2
    for set_time in [the_set_time, nothing]
        for i in eachindex(headers)
            try
                idum, time, nw, nh = EFIT.parse_gfile_header(headers[i]; set_time=set_time)
                if !isnothing(set_time)
                    @test time == the_set_time
                else
                    @test time == times[i]
                end
            catch
                if expect_exception[i]
                    println("Got the expected exception")
                else
                    println("Got the unexpected exception")
                end
                @test expect_exception[i]
            end
        end
    end
end

@testset "imas" begin
    dd = IMASdd.dd()
    eqt = resize!(dd.equilibrium.time_slice, 1)[1]
    EFIT.geqdsk2imas!(g, eqt; wall=dd.wall, add_derived=true)
    @test length(eqt.profiles_2d[1].grid.dim1) > 1
end
