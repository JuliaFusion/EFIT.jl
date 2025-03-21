"""
    geqdsk2imas!(
        gs::Vector{GEQDSKFile},
        dd::IMASdd.dd;
        geqdsk_cocos::Int=0,
        dd_cocos::Int=11,
        cocos_clockwise_phi::Bool=false,
        add_derived::Bool=false,
    )

Utility for writing equilibrium data from GEQDSKFiles to IMAS schema.

 * gs: Vector of GEQDSKFile instances
 * dd: Top level IMAS data dictionary
 * geqdsk_cocos: coordinate convention indentifier of the incoming GEQDSK set. 0 to auto detect.
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
 * cocos_clockwise_phi: hint for deciding COCOS, because GEQDSK probably doesn't constrain this.
 * add_derived: switch to use functions within EFIT.jl to extend equilibrium data
"""
function geqdsk2imas!(
    gs::Vector{GEQDSKFile},
    dd::IMASdd.dd;
    geqdsk_cocos::Int=0,
    dd_cocos::Int=11,
    cocos_clockwise_phi::Bool=false,
    add_derived::Bool=false,
)
    # Should work if dd is empty
    if ismissing(dd.dataset_description.data_entry, :pulse)
        dd.dataset_description.data_entry.pulse = parse(Int, split(split(gs[1].file, ".")[end-1], "g")[2])
    end
    geqdsk2imas!(
        gs, dd.equilibrium; wall=dd.wall,
        geqdsk_cocos=geqdsk_cocos, dd_cocos=dd_cocos, cocos_clockwise_phi=cocos_clockwise_phi,
        add_derived=add_derived,
    )
end

"""
    geqdsk2imas!(
        gs::Vector{GEQDSKFile},
        eq::IMASdd.equilibrium;
        wall=nothing,
        geqdsk_cocos::Int=0,
        dd_cocos::Int=11,
        cocos_clockwise_phi::Bool=false,
        add_derived::Bool=false,
    )

Utility for writing equilibrium data from GEQDSKFile instances to IMAS equilibrium IDS

 * gs: Vector of GEQDSKFile instances
 * eq: Reference to the equilibrium IDS
 * wall: Reference to the wall IDS
 * geqdsk_cocos: coordinate convention indentifier of the incoming GEQDSK set. 0 to auto detect.
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
 * cocos_clockwise_phi: hint for deciding COCOS, because GEQDSK probably doesn't constrain this.
 * add_derived: switch to use functions within EFIT.jl to extend equilibrium data
"""
function geqdsk2imas!(
    gs::Vector{GEQDSKFile},
    eq::IMASdd.equilibrium;
    wall=nothing,
    geqdsk_cocos::Int=0,
    dd_cocos::Int=11,
    cocos_clockwise_phi::Bool=false,
    add_derived::Bool=false,
)
    dd = IMASdd.top_dd(eq)
    if ismissing(dd.dataset_description.data_entry, :pulse)
        dd.dataset_description.data_entry.pulse = parse(Int, split(split(gs[1].file, ".")[end-1], "g")[2])
    end
    if wall !== nothing
        geqdsk2wall!(gs[1], wall; geqdsk_cocos=geqdsk_cocos, dd_cocos=dd_cocos, cocos_clockwise_phi=cocos_clockwise_phi)
    end
    nt = length(gs)
    if length(eq.time_slice) == 0  # Blank is okay, but if populated, it had better be compatible
        resize!(eq.time_slice, nt)
        eq.time = [g.time for g in gs]
        for i in 1:nt
            eq.time_slice[i].time = eq.time[i]
        end
    end

    # Write top level data
    eq.time = [g.time for g in gs]
    bb = zeros(length(gs))
    for i in 1:length(gs)
        if geqdsk_cocos == 0
            g = gs[i]
            geqdsk_cocos_ = CoordinateConventions.identify_cocos(
                sign(g.bcentr), sign(g.current), sign(g.qpsi[1]), sign(g.psi[end] - g.psi[1]), cocos_clockwise_phi,
            )[1]
        else
            geqdsk_cocos_ = geqdsk_cocos
        end
        tc = CoordinateConventions.transform_cocos(geqdsk_cocos_, dd_cocos)
        bb[i] = g.bcentr .* tc["B"]
    end
    eq.vacuum_toroidal_field.b0 = bb
    eq.vacuum_toroidal_field.r0 = gs[1].rcentr

    for g in gs
        geqdsk2imas!(
            g, eq;
            geqdsk_cocos=geqdsk_cocos, dd_cocos=dd_cocos, cocos_clockwise_phi=cocos_clockwise_phi,
            add_derived=add_derived,
        )
    end
end

"""
    geqdsk2imas!(
        g::GEQDSKFile,
        dd::IMASdd.dd;
        geqdsk_cocos::Int=0,
        dd_cocos::Int=11,
        cocos_clockwise_phi::Bool=false,
        add_derived::Bool=false,
    )

Utility for writing equilibrium data from a GEQDSKFile to IMAS DD.

 * g: GEQDSKFile instance
 * dd: IMAS data dictionary
 * geqdsk_cocos: coordinate convention indentifier of the incoming GEQDSK set. 0 to auto detect.
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
 * cocos_clockwise_phi: hint for deciding COCOS, because GEQDSK probably doesn't constrain this.
 * add_derived: switch to use functions within EFIT.jl to extend equilibrium data
"""
function geqdsk2imas!(
    g::GEQDSKFile,
    dd::IMASdd.dd;
    geqdsk_cocos::Int=0,
    dd_cocos::Int=11,
    cocos_clockwise_phi::Bool=false,
    add_derived::Bool=false,
)
    if ismissing(dd.dataset_description.data_entry, :pulse)
        dd.dataset_description.data_entry.pulse = parse(Int, split(split(g.file, ".")[end-1], "g")[2])
    end
    geqdsk2imas!(
        g, dd.equilibrium;
        dd.wall, geqdsk_cocos, dd_cocos, cocos_clockwise_phi, add_derived,
    )
end

"""
    geqdsk2imas!(
        g::GEQDSKFile,
        eq::IMASdd.equilibrium;
        wall=nothing,
        geqdsk_cocos::Int=0,
        dd_cocos::Int=11,
        cocos_clockwise_phi::Bool=false,
        add_derived::Bool=false,
    )

Utility for writing equilibrium data from a GEQDSKFile to IMAS equilibrium IDS

 * g: GEQDSKFile instance
 * eq: Reference to equilibrium IDS
 * wall: Optional reference to wall IDS
 * geqdsk_cocos: coordinate convention indentifier of the incoming GEQDSK set. 0 to auto detect.
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
 * cocos_clockwise_phi: hint for deciding COCOS, because GEQDSK probably doesn't constrain this.
 * add_derived: switch to use functions within EFIT.jl to extend equilibrium data
"""
function geqdsk2imas!(
    g::GEQDSKFile,
    eq::IMASdd.equilibrium;
    wall=nothing,
    geqdsk_cocos::Int=0,
    dd_cocos::Int=11,
    cocos_clockwise_phi::Bool=false,
    add_derived::Bool=false,
)
    if geqdsk_cocos == 0
        geqdsk_cocos = CoordinateConventions.identify_cocos(
            sign(g.bcentr), sign(g.current), sign(g.qpsi[1]), sign(g.psi[end] - g.psi[1]), cocos_clockwise_phi,
        )[1]
        @debug "  identified COCOS=$geqdsk_cocos"
    end
    tc = CoordinateConventions.transform_cocos(geqdsk_cocos, dd_cocos)

    dd = IMASdd.top_dd(eq)
    original_global_time = dd.global_time
    try
        dd.global_time = g.time
        resize!(eq.time_slice)
        geqdsk2imas!(
            g, eq.time_slice[];
            wall, geqdsk_cocos, dd_cocos, cocos_clockwise_phi, add_derived,
        )
    finally
        dd.global_time = original_global_time
    end
end

"""
    geqdsk2imas!(
        g::GEQDSKFile,
        eqt::IMASdd.equilibrium__time_slice{Float64};
        wall=nothing,
        geqdsk_cocos::Int=0,
        dd_cocos::Int=11,
        cocos_clockwise_phi::Bool=false,
        add_derived::Bool=false,
    )

Writes equilibrium data from a GEQDSK file to IMAS equilibrium IDS in a specific
time slice. Can optionally include writing wall data to the wall IDS.

 * g: GEQDSKFile instance
 * eqt: A specific time slice within dd.equilibrium.time_slice
 * wall: Optional reference to wall IDS
 * geqdsk_cocos: coordinate convention indentifier of the incoming GEQDSK set. 0 to auto detect.
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
 * cocos_clockwise_phi: hint for deciding COCOS, because GEQDSK probably doesn't constrain this.
 * add_derived: switch to use functions within EFIT.jl to extend equilibrium data
"""
function geqdsk2imas!(
    g::GEQDSKFile,
    eqt::IMASdd.equilibrium__time_slice{Float64};
    wall=nothing,
    geqdsk_cocos::Int=0,
    dd_cocos::Int=11,
    cocos_clockwise_phi::Bool=false,
    add_derived::Bool=false,
)
    if geqdsk_cocos == 0
        geqdsk_cocos = CoordinateConventions.identify_cocos(
            sign(g.bcentr), sign(g.current), sign(g.qpsi[1]), sign(g.psi[end] - g.psi[1]), cocos_clockwise_phi,
        )[1]
        @debug "  identified COCOS=$geqdsk_cocos"
    end
    tc = CoordinateConventions.transform_cocos(geqdsk_cocos, dd_cocos)

    dd = IMASdd.top_dd(eqt)
    original_global_time = dd.global_time
    try
        dd.global_time = g.time
        IMASdd.@ddtime(dd.equilibrium.vacuum_toroidal_field.b0 = g.bcentr .* tc["B"])
        dd.equilibrium.vacuum_toroidal_field.r0 = g.rcentr
    finally
        dd.global_time = original_global_time
    end

    # Global and boundary
    gq = eqt.global_quantities
    gq.magnetic_axis.r = g.rmaxis * tc["R"]
    gq.magnetic_axis.z = g.zmaxis * tc["Z"]
    gq.ip = g.current * tc["I"]
    gq.psi_axis = g.simag * tc["PSI"]
    gq.psi_boundary = g.sibry * tc["PSI"]
    eqt.boundary.geometric_axis.r = g.rcentr * tc["R"]
    eqt.boundary.geometric_axis.z = g.zmid * tc["Z"]
    eqt.boundary.outline.r = g.rbbbs .* tc["R"]
    eqt.boundary.outline.z = g.zbbbs .* tc["Z"]

    # 1D profiles
    p1 = eqt.profiles_1d
    p1.psi = g.psi .* tc["PSI"]
    p1.q = g.qpsi .* tc["Q"]
    p1.pressure = g.pres .* tc["P"]
    p1.dpressure_dpsi = g.pprime .* tc["PPRIME"]
    p1.f = g.fpol .* tc["F"]
    p1.f_df_dpsi = g.ffprim .* tc["F_FPRIME"]
    if hasproperty(g, :rhovn)
        p1.rho_tor_norm = g.rhovn
    end

    # 2D flux map
    p2 = resize!(eqt.profiles_2d, 1)[1]
    p2.grid_type.index = 1
    p2.grid_type.name = "R-Z grid for flux map"
    p2.grid_type.description = (
        "A rectangular grid of points in R,Z on which poloidal " *
        "magnetic flux psi is defined. The grid's dim1 is R, dim2 is Z."
    )
    p2.grid.dim1 = collect(g.r) .* tc["R"]
    p2.grid.dim2 = collect(g.z) .* tc["Z"]
    p2.psi = g.psirz .* tc["PSI"]

    # Wall
    if wall !== nothing
        geqdsk2wall!(g, wall; geqdsk_cocos=geqdsk_cocos, dd_cocos=dd_cocos)
    end

    # Derived
    if add_derived
        derived_g2imas!(g, eqt; geqdsk_cocos=geqdsk_cocos, dd_cocos=dd_cocos)
    end
end

"""
    derived_g2imas!(
        g::GEQDSKFile,
        eqt::IMASdd.equilibrium__time_slice;
        geqdsk_cocos::Int=0,
        dd_cocos::Int=11,
        cocos_clockwise_phi::Bool=false,
    )

Does simple calculations related to the flux map and stores results in IMAS

 * g: GEQDSKFile instance
 * eqt: A specific time slice within dd.equilibrium.time_slice
 * geqdsk_cocos: coordinate convention indentifier of the incoming GEQDSK set. 0 to auto detect.
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
 * cocos_clockwise_phi: hint for deciding COCOS, because GEQDSK probably doesn't constrain this.
"""
function derived_g2imas!(
    g::GEQDSKFile,
    eqt::IMASdd.equilibrium__time_slice;
    geqdsk_cocos::Int=0,
    dd_cocos::Int=11,
    cocos_clockwise_phi::Bool=false,
)
    if geqdsk_cocos == 0
        geqdsk_cocos = CoordinateConventions.identify_cocos(
            sign(g.bcentr), sign(g.current), sign(g.qpsi[1]), sign(g.psi[end] - g.psi[1]), cocos_clockwise_phi
        )[1]
        @debug "  identified COCOS=$geqdsk_cocos"
    end
    tc = CoordinateConventions.transform_cocos(geqdsk_cocos, dd_cocos)

    # X-points
    xrs, xzs, xpsins, xseps = x_points(g; within_limiter_only=false)
    xrs .*= tc["R"]
    xzs .*= tc["Z"]
    if length(xrs) > 0
        bx = eqt.boundary.x_point
        resize!(bx, length(xrs))
        for i ∈ eachindex(xrs)
            bx[i].r = xrs[i]
            bx[i].z = xzs[i]
        end
        nprim = sum(xseps .== 1)
        if nprim > 0
            bsx = eqt.boundary_separatrix.x_point
            resize!(bsx, nprim)
            xrprim = xrs[xseps.==1]
            xzprim = xzs[xseps.==1]
            for i ∈ nprim
                bsx[i].r = xrprim[i]
                bsx[i].z = xzprim[i]
            end
        end
        nsec = sum(xseps .== 2)
        if nsec > 0
            bssx = eqt.boundary_secondary_separatrix.x_point
            resize!(bssx, nsec)
            xrsec = xrs[xseps.==2]
            xzsec = xzs[xseps.==2]
            for i ∈ nsec
                bssx[i].r = xrsec[i]
                bssx[i].z = xzsec[i]
            end
        end
    end
end

"""
    geqdsk2wall!(
        g::GEQDSKFile,
        wall::IMASdd.wall;
        geqdsk_cocos::Int=0,
        dd_cocos::Int=11,
        cocos_clockwise_phi::Bool=false,
    )

Writes wall data from GEQDSK to the wall IDS in IMAS.

 * g: GEQDSKFile instance
 * wall: Reference to the wall IDS in IMAS DD
 * geqdsk_cocos: coordinate convention indentifier of the incoming GEQDSK set. 0 to auto detect.
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
 * cocos_clockwise_phi: hint for deciding COCOS, because GEQDSK probably doesn't constrain this.
"""
function geqdsk2wall!(
    g::GEQDSKFile,
    wall::IMASdd.wall;
    geqdsk_cocos::Int=0,
    dd_cocos::Int=11,
    cocos_clockwise_phi::Bool=false,
)
    if geqdsk_cocos == 0
        geqdsk_cocos = CoordinateConventions.identify_cocos(
            sign(g.bcentr), sign(g.current), sign(g.qpsi[1]), sign(g.psi[end] - g.psi[1]), cocos_clockwise_phi,
        )[1]
        @debug "  identified COCOS=$geqdsk_cocos"
    end
    tc = CoordinateConventions.transform_cocos(geqdsk_cocos, dd_cocos)
    resize!(wall.description_2d, 1)
    limiter = wall.description_2d[1].limiter
    limiter.type.name = "first wall"
    limiter.type.index = 0
    limiter.type.description = "first wall"
    resize!(limiter.unit, 1)
    limiter.unit[1].outline.r = g.rlim .* tc["R"]
    limiter.unit[1].outline.z = g.zlim .* tc["Z"]
end

"""
    imas2geqdsk(
        dd::IMASdd.dd;
        geqdsk_cocos::Int=1,
        dd_cocos::Int=11,
    )::Vector{GEQDSKFile}

Utility for outputting IMAS equilibrium data to a vector of GEQDSKFile instances.

 * dd: IMAS data dictionary instance containing equilibrium data
 * geqdsk_cocos: coordinate convention indentifier of the GEQDSKs to be written
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
"""
function imas2geqdsk(
    dd::IMASdd.dd;
    geqdsk_cocos::Int=1,
    dd_cocos::Int=11,
)::Vector{GEQDSKFile}
    nt = length(dd.equilibrium.time_slice)
    return [imas2geqdsk(dd, time, geqdsk_cocos=geqdsk_cocos, dd_cocos=dd_cocos) for time in dd.equilibrium.time]
end

"""
    imas2geqdsk(
        dd::IMASdd.dd,
        time::Float64;
        geqdsk_cocos::Int=1,
        dd_cocos::Int=11,
    )::GEQDSKFile

Utility for outputting IMAS equilibrium data from a single time slice to a GEQDSKFile

 * dd: IMAS data dictionary instance containing equilibrium data
 * time: Time of the slice to output in s. Lookup using IMAS global time utilities.
 * geqdsk_cocos: coordinate convention indentifier of the GEQDSKs to be written
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
"""
function imas2geqdsk(
    dd::IMASdd.dd,
    time::Float64;
    geqdsk_cocos::Int=1,
    dd_cocos::Int=11,
)::GEQDSKFile
    tc = CoordinateConventions.transform_cocos(geqdsk_cocos, dd_cocos)

    if !ismissing(dd.dataset_description.data_entry, :pulse)
        shot = dd.dataset_description.data_entry.pulse
    else
        shot = 0
    end

    original_global_time = dd.global_time
    try
        dd.global_time = time
        eqt = dd.equilibrium.time_slice[]
        limiter = dd.wall.description_2d[1].limiter
        rlim = limiter.unit[1].outline.r ./ tc["R"]
        zlim = limiter.unit[1].outline.z ./ tc["Z"]
        limitr = length(rlim)

        # Global and boundary
        gq = eqt.global_quantities
        rmaxis = gq.magnetic_axis.r / tc["R"]
        zmaxis = gq.magnetic_axis.z / tc["Z"]
        current = gq.ip / tc["I"]
        simag = gq.psi_axis / tc["PSI"]
        sibry = gq.psi_boundary / tc["PSI"]
        rcentr = eqt.boundary.geometric_axis.r / tc["R"]
        zmid = eqt.boundary.geometric_axis.z / tc["Z"]
        rbbbs = eqt.boundary.outline.r ./ tc["R"]
        zbbbs = eqt.boundary.outline.z ./ tc["Z"]
        nbbbs = length(rbbbs)

        # 2D flux map
        p2 = findfirst(:rectangular, eqt.profiles_2d)
        r = p2.grid.dim1 ./ tc["R"]
        z = p2.grid.dim2 ./ tc["Z"]
        psirz = p2.psi ./ tc["PSI"]

        # 1D profiles
        function interpolate_1d_profile(ori_1D::AbstractVector{<:Real}, N::Int)
            if length(ori_1D) == N
                return ori_1D
            else
                itp=IMASdd.interp1d(p1.psi_norm, ori_1D, :cubic)
                return itp.(range(0,1,N))
            end
        end

        p1 = eqt.profiles_1d
        nw = length(r)
        psi = interpolate_1d_profile(p1.psi ./ tc["PSI"], nw)
        qpsi = interpolate_1d_profile(p1.q ./ tc["Q"], nw)
        pres = interpolate_1d_profile(p1.pressure ./ tc["P"], nw)
        pprime = interpolate_1d_profile(p1.dpressure_dpsi ./ tc["PPRIME"], nw)
        fpol = interpolate_1d_profile(p1.f ./ tc["F"], nw)
        ffprim = interpolate_1d_profile(p1.f_df_dpsi ./ tc["F_FPRIME"], nw)
        rhovn = interpolate_1d_profile(p1.rho_tor_norm, nw)

        bcentr = IMASdd.@ddtime(dd.equilibrium.vacuum_toroidal_field.b0) ./ tc["B"]
        time = eqt.time
        @debug "eqt.time = $(eqt.time), time out = $time"
        rleft = minimum(r)
        rdim = maximum(r) - rleft
        zdim = maximum(z) - minimum(z)
        nw = length(r)
        nh = length(z)
        itime = Int(round(time * 1000))
        gfile = "g" * lpad(shot, 6, '0') * "." * lpad(itime, 5, '0')
        @debug "shot $shot, gfile = $gfile"

        r = range(rleft, rleft + rdim, length=nw)
        z = range(zmid - 0.5*zdim, zmid + 0.5*zdim, length=nh)
        psi = range(simag, sibry, length=nw)

        return GEQDSKFile(
            gfile,time,nw,nh,r,z,rdim,zdim,rleft,zmid,nbbbs,rbbbs,zbbbs,limitr,rlim,zlim,
            rcentr,bcentr,rmaxis,zmaxis,simag,sibry,psi,current,fpol,pres,ffprim,pprime,
            qpsi,psirz,rhovn,
        )
    finally
        dd.global_time = original_global_time
    end
end