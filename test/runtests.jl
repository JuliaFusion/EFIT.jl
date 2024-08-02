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

function gcompare(g1::GEQDSKFile, g2::GEQDSKFile; note=nothing, verbose::Bool=false)
    fields = split("file,time,nw,nh,r,z,rdim,zdim,rleft,zmid,nbbbs,rbbbs,zbbbs,limitr,rlim,zlim,"*
    "rcentr,bcentr,rmaxis,zmaxis,simag,sibry,psi,current,fpol,pres,ffprim,pprime,"*
    "qpsi,psirz,rhovn", ",")
    mismatched_fields = []
    if note !== nothing
        println(note)
    end
    for field in fields
        if verbose
            print(" $field: ")
        end
        f1 = getproperty(g1, Symbol(field))
        if verbose
            if (length(f1) > 3) && (typeof(f1) !== String)
                print(" [$(length(f1)) elements], ")
            else
                print(" $f1, ")
            end
        end
        f2 = getproperty(g2, Symbol(field))
        if verbose
            if (length(f1) > 3) && (typeof(f1) !== String)
                println("[$(length(f2)) elements]")
            else
                println(f2)
            end
        end
        if (typeof(f1) == String) || (length(f1) !== length(f2))
            if f1 !== f2
                mismatched_fields = [mismatched_fields; field]
            end
            @test f1 == f2
        else
            if !isapprox(f1, f2, atol=1e-5)
                mismatched_fields = [mismatched_fields; field]
            end
            @test f1 ≈ f2 atol=1e-5
        end
    end
    println("mismatched fields for this case = $mismatched_fields")
    return mismatched_fields
end

@testset "imas" begin
    gs = [g2, g]
    # Add one geqdsk
    println("test file 1 info: filename = $(gs[1].file), time = $(gs[1].time)")
    dd = IMASdd.dd()
    eqt = resize!(dd.equilibrium.time_slice, 1)[1]
    EFIT.geqdsk2imas!(gs[1], eqt; wall=dd.wall, add_derived=true)
    @test length(eqt.profiles_2d[1].grid.dim1) > 1
    @test eqt.time == gs[1].time
    # Add another geqdsk to another index
    println("test file 2 info: filename = $(gs[2].file), time = $(gs[2].time)")
    idx = 2
    EFIT.geqdsk2imas!(gs[2], dd.equilibrium, idx, add_derived=true)
    eqt2 = dd.equilibrium.time_slice[idx]
    @test length(eqt2.profiles_2d[1].grid.dim1) > 1
    @test abs(dd.equilibrium.vacuum_toroidal_field.b0[idx]) > 0.0
    @test eqt2.time == gs[2].time
    # Add multiple geqdsks
    newdd = IMASdd.dd()
    EFIT.geqdsk2imas!(gs, newdd; add_derived=true)
    @test length(newdd.equilibrium.time_slice) == length(gs)
    shot = parse(Int, split(split(gs[1].file, "g")[end], ".")[1])
    @test newdd.dataset_description.data_entry.pulse == shot
    println("newdd pulse is set to $(newdd.dataset_description.data_entry.pulse)")
    for slice in 1:length(gs)
        eqt1 = dd.equilibrium.time_slice[1]
        eqt2 = newdd.equilibrium.time_slice[1]
        for field in [:ip, :psi_axis, :psi_boundary]
            @test getproperty(eqt1.global_quantities, field) == getproperty(eqt2.global_quantities, field)
        end
        @test newdd.equilibrium.time[slice] == gs[slice].time
    end

    # Reverse the process
    gs2 = EFIT.imas2geqdsk(newdd)
    gg1 = EFIT.imas2geqdsk(newdd, 1)
    gg2 = EFIT.imas2geqdsk(newdd, 2)
    # There are some artifacts of using a set of random g-files for testing within a single
    # dd instead of a coherent set of slices from the same shot.
    # Only the first shot is used.
    # The contents of the file have different time units, making it hard to get the right
    # filename from the time. So, we'll reset the filename of the second file, which is also
    # the one with weird time units.
    gg2.file = gs2[2].file= gs[2].file
    # Only the wall from the first file is used.
    gg2.rlim = gs2[2].rlim = gs[2].rlim
    gg2.zlim = gs2[2].zlim = gs[2].zlim
    gg2.limitr = gs2[2].limitr = gs[2].limitr
    gcompare(gg1, gs2[1], note="read single g vs read all g from imas (slice 1)")
    gcompare(gg2, gs2[2], note="read single g vs read all g from imas (slice 2)")
    gcompare(gg1, gs[1], note="read single g vs input 1")
    gcompare(gg2, gs[2], note="read single g vs input 2")

    # Make sure something weird about the second file doesn't break the wall situation
    dd3 = IMASdd.dd()
    EFIT.geqdsk2imas!(gs[2], dd3, 1)
    gg3 = EFIT.imas2geqdsk(dd3, 1)
    # as this is the sample with the funky time units in the header, we'll help with the filename
    gg3.file = gs[2].file
    gcompare(gg3, gs[2])
end
