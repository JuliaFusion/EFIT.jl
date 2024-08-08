import IMASdd
import CoordinateConventions

mutable struct GEQDSKFile
    file::String                    # Source file
    time::Float64                   # Time of the equilibrium reconstruction, in seconds
    nw::Int                         # Number of horizontal R grid points
    nh::Int                         # Number of vertical Z grid points
    r::AbstractRange{Float64}       # R grid points
    z::AbstractRange{Float64}       # Z grid points
    rdim::Float64                   # Horizontal dimension in meter of computational box
    zdim::Float64                   # Vertical dimension in meter of computational box
    rleft::Float64                  # Minimum R in meter of rectangular computational box
    zmid::Float64                   # Z of center of computational box in meter
    nbbbs::Int                      # Number of boundary points
    rbbbs::Vector{Float64}          # R of boundary points in meter
    zbbbs::Vector{Float64}          # Z of boundary points in meter
    limitr::Int                     # Number of limiter points
    rlim::Vector{Float64}           # R of surrounding limiter contour in meter
    zlim::Vector{Float64}           # Z of surrounding limiter contour in meter
    rcentr::Float64                 # R in meter of vacuum toroidal magnetic field BCENTR
    bcentr::Float64                 # Vacuum toroidal magnetic field in Tesla at RCENTR
    rmaxis::Float64                 # R of magnetic axis in meter
    zmaxis::Float64                 # Z of magnetic axis in meter
    simag::Float64                  # Poloidal flux at magnetic axis in Weber/rad
    sibry::Float64                  # Poloidal flux at the plasma boundary in Weber/rad
    psi::AbstractRange{Float64}     # Poloidal flux points
    current::Float64                # Plasma current in Ampere
    fpol::Vector{Float64}           # Poloidal current function in m-T, F = RBT on flux grid
    pres::Vector{Float64}           # Plasma pressure in nt / m2 on uniform flux grid
    ffprim::Vector{Float64}         # FF'(psi) in (mT)2 / (Weber/rad) on uniform fl ux grid
    pprime::Vector{Float64}         # P'(psi) in (nt/m2) / (Weber/rad) on uniform flux grid
    qpsi::Vector{Float64}           # q values on uniform flu x grid from axis to boundary
    psirz::Matrix{Float64}          # Poloidal flux in Weber/rad on the rectangular grid points
    rhovn::Vector{Float64}          # Square root of toroidal flux, normalized
end

function Base.show(io::IO, g::GEQDSKFile)
    print(io,"GEQDSKFile: \"",g.file,"\"")
end

Base.broadcastable(g::GEQDSKFile) = (g,)

function file_numbers(f)
    Channel(ctype=Float64) do c
        while ~eof(f)
            line = readline(f)
            pattern = r"[+-]?\d*[\.]?\d+(?:[Ee][+-]?\d+)?"
            for p in eachmatch(pattern,line)
                push!(c,parse(Float64,p.match))
            end
        end
    end
end

function read_array(t,n)
    data = zeros(n)
    for i=1:n
        try
            data[i] = take!(t)
        catch e
            if isa(e, InvalidStateException)
                # InvalidStateException when Channel is closed
            else
                rethrow(e)
            end
        end
    end
    return data
end

function read_array2d(t,n,m)
    data = zeros(n,m)
    for i=1:n, j=1:m
        try
            data[i,j] = take!(t)
        catch e
            if isa(e, InvalidStateException)
                # InvalidStateException when Channel is closed
            else
                rethrow(e)
            end
        end
    end
    return data'
end

function parse_gfile_header(headerline::String; set_time=nothing)
    s = split(headerline)

    if !isnothing(set_time)
        time = set_time
    else
        time = Meta.parse(s[end-3])
        time_warning = "EFIT.jl could not convert time from ms (assumed units) to s"
        omfit_advice = "Options for proceeding include:\n" *
                       "i)    Editing the header in your gEQDSK file.\n" *
                       "ii)   Trying to use OMFIT: python: `OMFITgeqdsk(file).to_omas().save('ods.jl')`\n" *
                       "iii)  Set the time manually with the `set_time` keyword argument.\n"
        try
            time /= 1000.0  # Assume it was written in ms and convert to s
        catch e
            if time isa String
                rg2 = r"([[:digit:]]+|(?:[[:punct:]]|[[:blank:]])+)"
                try
                    time = parse(Int, collect(eachmatch(rg2, time))[1][1]) / 1000.0
                catch e
                    error(
                        time_warning,
                        ", even after trying really hard with regex and everything.\nTime is ",
                        time,
                        ". \n",
                        omfit_advice
                    )
                end
            else
                error(time_warning, ".\nTime is parsed as ", time," of type ", typeof(time), ". ", omfit_advice)
            end
        end
    end

    idum = Meta.parse(s[end-2])
    nw = Meta.parse(s[end-1])
    nh = Meta.parse(s[end])

    return idum, time, nw, nh
end

function readg(gfile; set_time=nothing)

    isfile(gfile) || error("$(gfile) does not exist")

    f = open(gfile)

    desc = readline(f)
    idum, time, nw, nh = parse_gfile_header(desc; set_time)

    token = file_numbers(f)

    rdim   = take!(token)
    zdim   = take!(token)
    rcentr = take!(token)
    rleft  = take!(token)
    zmid   = take!(token)

    rmaxis = take!(token)
    zmaxis = take!(token)
    simag  = take!(token)
    sibry  = take!(token)
    bcentr = take!(token)

    current = take!(token)
    simag  = take!(token)
    xdum    = take!(token)
    rmaxis  = take!(token)
    xdum    = take!(token)

    zmaxis = take!(token)
    xdum   = take!(token)
    sibry = take!(token)
    xdum   = take!(token)
    xdum   = take!(token)

    fpol = read_array(token,nw)
    pres = read_array(token,nw)
    ffprim = read_array(token,nw)
    pprime = read_array(token,nw)
    psirz = read_array2d(token,nh,nw)
    qpsi = read_array(token,nw)

    nbbbs = Int(take!(token))
    limitr = Int(take!(token))

    if nbbbs > 0
        rbbbs = zeros(nbbbs)
        zbbbs = zeros(nbbbs)
        for i=1:nbbbs
            rbbbs[i] = take!(token)
            zbbbs[i] = take!(token)
        end
    else
        rbbbs = [0.0]
        zbbbs = [0.0]
    end

    if limitr > 0
        rlim = zeros(limitr)
        zlim = zeros(limitr)
        for i=1:limitr
            rlim[i] = take!(token)
            zlim[i] = take!(token)
        end
    else
        rlim = [0.0]
        zlim = [0.0]
    end

    local rhovn
    try
        for i ∈ 1:3
            if !eof(f)
                xdum = take!(token)
            end
        end
        if !eof(f)
            rhovn = read_array(token,nw)
        else
            rhovn = zeros(nw)
        end
        if !eof(f)
            xdum = take!(token)
        end
    catch e
        if isa(e, InvalidStateException)
            # InvalidStateException when Channel is closed
            rhovn = zeros(nw)
        else
            rethrow(e)
        end
    end

    close(f)
    close(token)

    r = range(rleft, rleft + rdim, length=nw)
    z = range(zmid - 0.5*zdim, zmid + 0.5*zdim, length=nh)
    psi = range(simag, sibry, length=nw)

    g = GEQDSKFile(
        splitpath(gfile)[end],time,nw,nh,r,z,rdim,zdim,rleft,zmid,nbbbs,rbbbs,zbbbs,
        limitr,rlim,zlim,rcentr,bcentr,rmaxis,zmaxis,simag,sibry,psi,current,fpol,pres,
        ffprim,pprime,qpsi,psirz,rhovn
    )

    return g
end

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
    if wall !== nothing
        geqdsk2wall!(gs[1], wall; geqdsk_cocos=geqdsk_cocos, dd_cocos=dd_cocos, cocos_clockwise_phi=cocos_clockwise_phi)
    end
    nt = length(gs)
    if length(eq.time_slice) < nt
        resize!(eq.time_slice, nt)
    end
    for it in 1:nt
        g = gs[it]
        geqdsk2imas!(
            g, eq, it;
            geqdsk_cocos=geqdsk_cocos, dd_cocos=dd_cocos, cocos_clockwise_phi=cocos_clockwise_phi,
            add_derived=add_derived,
        )
    end
end

"""
    geqdsk2imas!(
        g::GEQDSKFile,
        dd::IMASdd.dd,
        time_index::Int;
        geqdsk_cocos::Int=0,
        dd_cocos::Int=11,
        cocos_clockwise_phi::Bool=false,
        add_derived::Bool=false,
    )

Utility for writing equilibrium data from a GEQDSKFile to IMAS DD.

 * g: GEQDSKFile instance
 * dd: IMAS data dictionary
 * time_index: index of the destination for equilibrium data within equilibrium.time_slice
 * geqdsk_cocos: coordinate convention indentifier of the incoming GEQDSK set. 0 to auto detect.
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
 * cocos_clockwise_phi: hint for deciding COCOS, because GEQDSK probably doesn't constrain this.
 * add_derived: switch to use functions within EFIT.jl to extend equilibrium data
"""
function geqdsk2imas!(
    g::GEQDSKFile,
    dd::IMASdd.dd,
    time_index::Int;
    geqdsk_cocos::Int=0,
    dd_cocos::Int=11,
    cocos_clockwise_phi::Bool=false,
    add_derived::Bool=false,
)
    if ismissing(dd.dataset_description.data_entry, :pulse)
        dd.dataset_description.data_entry.pulse = parse(Int, split(split(g.file, ".")[end-1], "g")[2])
    end
    geqdsk2imas!(
        g, dd.equilibrium, time_index;
        wall=dd.wall, geqdsk_cocos=geqdsk_cocos, dd_cocos=dd_cocos, cocos_clockwise_phi=cocos_clockwise_phi,
        add_derived=add_derived,
    )
end

"""
    geqdsk2imas!(
        g::GEQDSKFile,
        eq::IMASdd.equilibrium,
        time_index::Int;
        wall=nothing,
        geqdsk_cocos::Int=0,
        dd_cocos::Int=11,
        cocos_clockwise_phi::Bool=false,
        add_derived::Bool=false,
    )

Utility for writing equilibrium data from a GEQDSKFile to IMAS equilibrium IDS

 * g: GEQDSKFile instance
 * eq: Reference to equilibrium IDS
 * time_index: index of the destination for equilibrium data within equilibrium.time_slice
 * wall: Optional reference to wall IDS
 * geqdsk_cocos: coordinate convention indentifier of the incoming GEQDSK set. 0 to auto detect.
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
 * cocos_clockwise_phi: hint for deciding COCOS, because GEQDSK probably doesn't constrain this.
 * add_derived: switch to use functions within EFIT.jl to extend equilibrium data
"""
function geqdsk2imas!(
    g::GEQDSKFile,
    eq::IMASdd.equilibrium,
    time_index::Int;
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

    # Make sure top level stuff is defined
    if ismissing(eq, :time)
        eq.time = Array{Float64}(undef, time_index)
    end
    if time_index > length(eq.time)
        resize!(eq.time, time_index)
    end
    if ismissing(eq.vacuum_toroidal_field, :b0)
        eq.vacuum_toroidal_field.b0 = Array{Float64}(undef, time_index)
    end
    if time_index > length(eq.vacuum_toroidal_field.b0)
        resize!(eq.vacuum_toroidal_field.b0, time_index)
    end

    # Write top level data
    eq.time[time_index] = g.time
    eq.vacuum_toroidal_field.b0[time_index] = g.bcentr .* tc["B"]
    eq.vacuum_toroidal_field.r0 = g.rcentr

    # Handle time slice data
    if time_index > length(eq.time_slice)
        resize!(eq.time_slice, time_index)
    end
    eqt = eq.time_slice[time_index]
    geqdsk2imas!(
        g, eqt;
        wall=wall, geqdsk_cocos=geqdsk_cocos, dd_cocos=dd_cocos, cocos_clockwise_phi=cocos_clockwise_phi,
        add_derived=add_derived,
    )
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

    eqt.time = g.time

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
    return [imas2geqdsk(dd, time_index, geqdsk_cocos=geqdsk_cocos, dd_cocos=dd_cocos) for time_index in 1:nt]
end

"""
    imas2geqdsk(
        dd::IMASdd.dd,
        time_index::Int;
        geqdsk_cocos::Int=1,
        dd_cocos::Int=11,
    )::GEQDSKFile

Utility for outputting IMAS equilibrium data from a single time slice to a GEQDSKFile

 * dd: IMAS data dictionary instance containing equilibrium data
 * time_index: index of the time slice to output
 * geqdsk_cocos: coordinate convention indentifier of the GEQDSKs to be written
 * dd_cocos: coordinate convention identifier of the DD. Should normally be left at 11.
"""
function imas2geqdsk(
    dd::IMASdd.dd,
    time_index::Int;
    geqdsk_cocos::Int=1,
    dd_cocos::Int=11,
)::GEQDSKFile
    tc = CoordinateConventions.transform_cocos(geqdsk_cocos, dd_cocos)

    if !ismissing(dd.dataset_description.data_entry, :pulse)
        shot = dd.dataset_description.data_entry.pulse
    else
        shot = 0
    end

    eqt = dd.equilibrium.time_slice[time_index]
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

    # 1D profiles
    p1 = eqt.profiles_1d
    psi = p1.psi ./ tc["PSI"]
    qpsi = p1.q ./ tc["Q"]
    pres = p1.pressure ./ tc["P"]
    pprime = p1.dpressure_dpsi ./ tc["PPRIME"]
    fpol = p1.f ./ tc["F"]
    ffprim = p1.f_df_dpsi ./ tc["F_FPRIME"]
    rhovn = p1.rho_tor_norm

    # 2D flux map
    p2 = eqt.profiles_2d[1]
    r = p2.grid.dim1 ./ tc["R"]
    z = p2.grid.dim2 ./ tc["Z"]
    psirz = p2.psi ./ tc["PSI"]

    bcentr = dd.equilibrium.vacuum_toroidal_field.b0[time_index] ./ tc["B"]
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
end
