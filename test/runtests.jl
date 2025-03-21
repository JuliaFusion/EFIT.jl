using Test
using EFIT
import EFIT
import EFIT.IMASdd
import EFIT.CoordinateConventions

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

@testset "writeg" begin
    for test_gfile in [EFIT.test_gfile, EFIT.test_gfile2]

        ori_g = readg(EFIT.test_gfile)

        # Create and write a tmporary EQDSK file
        (tmp_filename, ~) = Base.mktemp()
        writeg(ori_g, tmp_filename; desc="\nHello, can you read me?\n\n")

        # Read the temporary EQDSK file
        new_g = readg(tmp_filename)

        # Match the "file" name for a fair comparison
        ori_g.file = "we got the same name for a fair comparison"
        new_g.file = "we got the same name for a fair comparison"

        # Check if new_g is equivalent to ori_g
        @test new_g == ori_g
        @test isequal(new_g, ori_g)

        # keywords test
        @test writeg(new_g, tmp_filename; desc="my description", shot="my shot 12345", time="1000ms")
        @test_throws Exception writeg(ori_g, tmp_filename; desc="This is a long description that will cause an error.")

        # Delete temporary file
        rm(tmp_filename, force=true)
    end
end

function gcompare(g1::GEQDSKFile, g2::GEQDSKFile; note=nothing, verbose::Bool=false)
    fields = fieldnames(GEQDSKFile)
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
    cocos_clockwise_phi = false  # Assumed
    gs = [g, g2]
    gcocos = [CoordinateConventions.identify_cocos(
        sign(gs[ig].bcentr), sign(gs[ig].current), sign(gs[ig].qpsi[1]), sign(gs[ig].psi[end] - gs[ig].psi[1]), cocos_clockwise_phi,
    )[1] for ig in 1:length(gs)]
    println("Test files are identified with COCOS = $gcocos")
    # Add one geqdsk
    println("test file 1 info: filename = $(gs[1].file), time = $(gs[1].time)")
    dd = IMASdd.dd()
    dd.equilibrium.time = [g.time for g in gs]

    resize!(dd.equilibrium.time_slice, length(gs))
    for (k, eqt) in pairs(dd.equilibrium.time_slice)
        eqt.time = gs[k].time
    end

    eqt = dd.equilibrium.time_slice[1]
    EFIT.geqdsk2imas!(gs[1], eqt; wall=dd.wall, cocos_clockwise_phi=cocos_clockwise_phi, add_derived=true)
    @test length(eqt.profiles_2d[1].grid.dim1) > 1
    @test eqt.time == gs[1].time
    # Add another geqdsk to another index
    println("test file 2 info: filename = $(gs[2].file), time = $(gs[2].time)")
    idx = 2
    dd.equilibrium.time_slice[2].time = gs[2].time
    EFIT.geqdsk2imas!(gs[2], dd.equilibrium, add_derived=true, cocos_clockwise_phi=cocos_clockwise_phi)
    eqt2 = dd.equilibrium.time_slice[idx]
    @test length(eqt2.profiles_2d[1].grid.dim1) > 1
    # @test abs(dd.equilibrium.vacuum_toroidal_field.b0[idx]) > 0.0
    @test eqt2.time == gs[2].time
    # Add multiple geqdsks
    newdd = IMASdd.dd()
    EFIT.geqdsk2imas!(gs, newdd; add_derived=true, cocos_clockwise_phi=cocos_clockwise_phi)
    @test length(newdd.equilibrium.time_slice) == length(gs)
    @test length(newdd.equilibrium.time) == length(gs)
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
        @test newdd.equilibrium.time_slice[slice].time == gs[slice].time
    end

    # Reverse the process
    # This part will not recover the original files if they have different cocos, since
    # the interface doesn't allow for per-file cocos specification. Per-file cocos
    # detection is possible on input, though.
    gs2 = EFIT.imas2geqdsk(newdd; geqdsk_cocos=gcocos[1])
    gg1 = EFIT.imas2geqdsk(newdd, gs[1].time; geqdsk_cocos=gcocos[1])
    gg2_orig = EFIT.imas2geqdsk(newdd, gs[2].time; geqdsk_cocos=gcocos[2])
    gg2_match = EFIT.imas2geqdsk(newdd, gs[2].time; geqdsk_cocos=gcocos[1])
    gg_out_orig = [gg1, gg2_orig]
    gg_out_match = [gg1, gg2_match]
    gcocos2 = [CoordinateConventions.identify_cocos(
        sign(gs2[ig].bcentr), sign(gs2[ig].current), sign(gs2[ig].qpsi[1]), sign(gs2[ig].psi[end] - gs2[ig].psi[1]), cocos_clockwise_phi,
    )[1] for ig in 1:length(gs2)]
    gcocos_out_orig = [CoordinateConventions.identify_cocos(
        sign(gg_out_orig[ig].bcentr), sign(gg_out_orig[ig].current), sign(gg_out_orig[ig].qpsi[1]), sign(gg_out_orig[ig].psi[end] - gg_out_orig[ig].psi[1]), cocos_clockwise_phi,
    )[1] for ig in 1:length(gg_out_orig)]
    gcocos_out_match = [CoordinateConventions.identify_cocos(
        sign(gg_out_match[ig].bcentr), sign(gg_out_match[ig].current), sign(gg_out_match[ig].qpsi[1]), sign(gg_out_match[ig].psi[end] - gg_out_match[ig].psi[1]), cocos_clockwise_phi,
    )[1] for ig in 1:length(gg_out_match)]
    for ig in 1:length(gs2)
        @test gcocos[ig] == gcocos_out_orig[ig]
        @test gcocos2[ig] == gcocos_out_match[ig]
    end
    # There are some artifacts of using a set of random g-files for testing within a single
    # dd instead of a coherent set of slices from the same shot.
    # Only the first shot is used.
    # The contents of the file have different time units, making it hard to get the right
    # filename from the time. So, we'll reset the filename of the second file, which is also
    # the one with weird time units.
    # One file has messed up time units
    gg1.file = gs2[1].file = gs[1].file
    # The second file won't have its original shot because the dd can only have one shot and it comes from the first one.
    gg2_match.file = gg2_orig.file = gs2[2].file= gs[2].file
    # Only the wall from the first file is used.
    gg2_orig.rlim = gg2_match.rlim = gs2[2].rlim = gs[2].rlim
    gg2_orig.zlim = gg2_match.zlim = gs2[2].zlim = gs[2].zlim
    gg2_orig.limitr = gg2_match.limitr = gs2[2].limitr = gs[2].limitr
    gcompare(gg1, gs2[1], note="read single g vs read all g from imas (slice 1)")
    gcompare(gg2_match, gs2[2], note="read single g vs read all g from imas (slice 2)")
    gcompare(gg1, gs[1], note="read single g vs input 1")
    gcompare(gg2_orig, gs[2], note="read single g vs input 2")

    # Make sure something weird about the second file doesn't break the wall situation
    dd3 = IMASdd.dd()
    EFIT.geqdsk2imas!(gs[2], dd3, cocos_clockwise_phi=cocos_clockwise_phi)
    gg3 = EFIT.imas2geqdsk(dd3, gs[2].time, geqdsk_cocos=gcocos[2])
    # as this is the sample with the funky time units in the header, we'll help with the filename
    gg3.file = gs[2].file
    gcompare(gg3, gs[2])
end

@testset "imas2geqdsk & geqdsk2imas!" begin
    ori_dd = IMASdd.json2imas("test_dd_eq.json")

    gg=EFIT.imas2geqdsk(ori_dd)

    new_dd = IMASdd.dd()
    EFIT.geqdsk2imas!(gg, new_dd)

    # compare wall
    @test new_dd.wall == ori_dd.wall
    @test new_dd.equilibrium.vacuum_toroidal_field == ori_dd.equilibrium.vacuum_toroidal_field

    ori_eqt = ori_dd.equilibrium.time_slice[]
    new_eqt = new_dd.equilibrium.time_slice[]

    # compare boundary
    @test new_eqt.boundary.geometric_axis == ori_eqt.boundary.geometric_axis
    @test new_eqt.boundary.outline == ori_eqt.boundary.outline

    # compare global_quantities
    ori_gq = new_eqt.global_quantities
    new_gq = new_eqt.global_quantities

    @test ori_gq.ip == new_gq.ip
    @test ori_gq.magnetic_axis == new_gq.magnetic_axis
    @test ori_gq.psi_axis == new_gq.psi_axis
    @test ori_gq.psi_boundary == new_gq.psi_boundary

    # compare profiles_1d
    ori_p1d = ori_eqt.profiles_1d
    new_p1d = new_eqt.profiles_1d

    target_fields = filter(x -> x ∉ IMASdd.private_fields, fieldnames(typeof(new_p1d)))
    for field in target_fields
        if !ismissing(new_p1d, field) && !isempty(new_p1d, field)
            ori1D_vec = getproperty(ori_p1d, field)
            new1D_vec = getproperty(new_p1d, field)

            itp=IMASdd.interp1d(ori_p1d.psi_norm, ori1D_vec, :cubic)
            @test isapprox(new1D_vec, itp.(range(0,1, length(new_p1d.psi))))
        end
    end

    # compare profiles_2d
    ori_p2d = findfirst(:rectangular, ori_eqt.profiles_2d)
    new_p2d = findfirst(:rectangular, new_eqt.profiles_2d)

    @test new_p2d.grid == ori_p2d.grid
    @test isapprox(new_p2d.psi, ori_p2d.psi)

end