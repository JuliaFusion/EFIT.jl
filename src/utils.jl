function triangularity(g::GEQDSKFile)
    Rmin, Rmax = extrema(g.rbbbs)
    Rgeo = (Rmin + Rmax)/2
    a = (Rmax - Rmin)/2
    Rupper = g.rbbbs[indmax(g.zbbbs)]
    Rlower = g.rbbbs[indmin(g.zbbbs)]
    delta_upper = (Rgeo - Rupper)/a
    delta_lower = (Rgeo - Rlower)/a

    return delta_lower, delta_upper
end

function ellipticity(g::GEQDSKFile)
    Rmin, Rmax = extrema(g.rbbbs)
    a = (Rmax - Rmin)/2
    Zmin, Zmax = extrema(g.zbbbs)

    return (Zmax - Zmin)/(2a)
end
