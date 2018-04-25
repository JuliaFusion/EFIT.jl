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
    return data
end

function readg(gfile)
    f = open(gfile)

    desc = readline(f)
    s = split(desc)

    time = parse(s[end-3])
    idum = parse(s[end-2])
    nw = parse(s[end-1])
    nh = parse(s[end])

    token = file_numbers(f)

    rdim   = take!(token)
    zdim   = take!(token)
    rcentr = take!(token)
    rleft  = take!(token)
    zmid   = take!(token)

    rmaxis = take!(token)
    zmaxis = take!(token)
    ssimag = take!(token)
    ssibry = take!(token)
    bcentr = take!(token)

    current = take!(token)
    ssimag  = take!(token)
    xdum    = take!(token)
    rmaxis  = take!(token)
    xdum    = take!(token)

    zmaxis = take!(token)
    xdum   = take!(token)
    ssibry = take!(token)
    xdum   = take!(token)
    xdum   = take!(token)

    fpol = read_array(token,nw)
    pres = read_array(token,nw)
    ffprim = read_array(token,nw)
    pprime = read_array(token,nw)
    psirz = read_array2d(token,nw,nh)
    qpsi = read_array(token,nw)

    nbdry = Int(take!(token))
    limitr = Int(take!(token))

    if nbdry > 0
        rbbbs = zeros(nbdry)
        zbbbs = zeros(nbdry)
        for i=1:nbdry
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

    close(f)
    close(token)

    r = linspace(rleft, rleft + rdim, nw)
    z = linspace(zmid - 0.5*zdim, zmid + 0.5*zdim, nh)

    g = Dict(:time=>time,:nw=>nw, :nh=>nh, :r=>r, :z=>z, :rdim=>rdim, :zdim=>zdim,
             :rcentr=>rcentr, :bcentr=>bcentr, :rleft=>rleft, :zmid=>zmid,
             :rmaxis=>rmaxis, :zmaxis=>zmaxis, :ssimag=>ssimag,:ssbry=>ssibry,
             :current=>current, :psirz=>psirz,:fpol=>fpol,:ffprim=>ffprim,:pprime=>pprime,
             :qpsi=>qpsi, :nbdry=>nbdry, :bdry=>hcat(rbbbs,zbbbs),
             :limitr=>limitr, :lim=>hcat(rlim,zlim))

    return g
end
