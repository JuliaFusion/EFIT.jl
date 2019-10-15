function triangularity(g::GEQDSKFile)
    Rmin, Rmax = extrema(filter(!iszero,g.rbbbs))
    Rgeo = (Rmin + Rmax)/2
    a = (Rmax - Rmin)/2
    Rupper = g.rbbbs[argmax(g.zbbbs)]
    Rlower = g.rbbbs[argmin(g.zbbbs)]
    delta_upper = (Rgeo - Rupper)/a
    delta_lower = (Rgeo - Rlower)/a

    return delta_lower, delta_upper
end

function ellipticity(g::GEQDSKFile)
    Rmin, Rmax = extrema(filter(!iszero,g.rbbbs))
    a = (Rmax - Rmin)/2
    Zmin, Zmax = extrema(g.zbbbs)

    return (Zmax - Zmin)/(2a)
end

function elongation(g::GEQDSKFile)
    ellipticity(g)
end

function major_radius(g::GEQDSKFile)
    0.5*(+(extrema(filter(!iszero,g.rbbbs))...))
end

function minor_radius(g::GEQDSKFile)
    -0.5*(-(extrema(filter(!iszero,g.rbbbs))...))
end

function aspect_ratio(g::GEQDSKFile)
    major_radius(g)/minor_radius(g)
end
