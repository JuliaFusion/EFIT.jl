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
        data[i] = take!(t)
    end
    return data
end

function read_array2d(t,n,m)
    data = zeros(n,m)
    for i=1:n, j=1:m
        data[i,j] = take!(t)
    end
    return data'
end

function readg(gfile)

    isfile(gfile) || error("$(gfile) does not exist")

    f = open(gfile)

    desc = readline(f)
    s = split(desc)

    time = Meta.parse(s[end-3])
    time /= 1000.0  # Assume it was written in ms and convert to s
    idum = Meta.parse(s[end-2])
    nw = Meta.parse(s[end-1])
    nh = Meta.parse(s[end])

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

    close(f)
    close(token)

    r = range(rleft, rleft + rdim, length=nw)
    z = range(zmid - 0.5*zdim, zmid + 0.5*zdim, length=nh)
    psi = range(simag, sibry, length=nw)

    g = GEQDSKFile(gfile,time,nw,nh,r,z,rdim,zdim,rleft,zmid,nbbbs,rbbbs,zbbbs,limitr,rlim,zlim,
                   rcentr,bcentr,rmaxis,zmaxis,simag,sibry,psi,current,fpol,pres,ffprim,pprime,
                   qpsi,psirz,rhovn)

    return g
end
