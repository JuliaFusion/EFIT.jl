using Interpolations: Interpolations
using PolygonOps: PolygonOps
using IMASdd: IMASdd

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

function triangularity(g::GEQDSKFile,ula::Symbol)
    δs = triangularity(g)
    if ula == :lower
        δ = δs[1]
    elseif ula == :upper
        δ = δs[2]
    elseif ula == :average
        δ = 0.5*sum(δs)
    else
        error("Unknown symbol ($ula): Only :lower, :upper, and :average accepted")
    end
    return δ
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

function elevation(g::GEQDSKFile)
    0.5*(+(extrema(g.zbbbs)...))
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

"""
    fluxinfo(r::Vector, z::Vector, psi::Matrix)

Calculates quantities derived from the flux map: R,Z, and polidal field components, and
discriminant for distinguishing minima and maxima from saddle points.

r: major radius in machine cylindrical coordinate system / m
z: height above midplane / m
psi: polidal magnetic flux / T m^2
"""
function fluxinfo(r::Vector, z::Vector, psi::Matrix)
    # Take some gradients
    dpsidr, dpsidz = IMASdd.gradient(r, z, psi)
    br = -dpsidz ./ r
    bz = dpsidr ./ r
    bpol2 = (br .^ 2) .+ (bz .^ 2)
    bpol = real.(sqrt.((br .^ 2) .+ (bz .^ 2)))
    d2psidr2, d2psidrdz = IMASdd.gradient(r, z, dpsidr)
    d2psidzdr, d2psidz2 = IMASdd.gradient(r, z, dpsidz)
    d = d2psidr2 .* d2psidz2 .- d2psidrdz .^ 2
    return (br, bz, bpol, d)
end

"""
    function find_local_minimum(a, neighborhood::Int64=3)

Finds the minimum within a neighborhood.
"""
function find_local_minimum(a, neighborhood::Int64=3)
    nx, ny = size(a)
    b = a * 0.0
    for i ∈ 1:nx
        imin = maximum([i - neighborhood, 1])
        imax = minimum([i + neighborhood, nx])
        for j ∈ 1:ny
            jmin = maximum([j - neighborhood, 1])
            jmax = minimum([j + neighborhood, ny])
            b[i, j] = minimum(a[imin:imax, jmin:jmax])
        end
    end
    return b
end

"""
    x_points(g::GEQDSKFile)

Given a flux map within a GEQDSKFile structure, find X-points. Return a list of
X-point R coordinates, Z coordinates, normalized psi values, and a list of flags
indicating to which separatrix (primary=1, secondary=2, other=4) they are
attached. Separatrix values other than 1, 2, or 4 indicate confusion.

Default settings hopefully work without adjustment, but if you need to fine tune:

neighborhood: number of mesh cells in either direction to include in considering
    a local neighborhood; this is for collapsing zones of low bpol into individual
    search regions centered on the local minimum of each region.
bpol_thresh_factor:
    One X-point clearly must be near the global minimum of bpol. But because of
    limited grid resolution, the local minimum of bpol near the next X-point
    will be larger than the global minimum. Some tolerance is needed to allow
    for searching near other low values of bpol.
psin_tolerance:
    Tolerance in psin for finding the X-point on the primary separatrix (should
    nominally have psin=1.0 exactly, by definition, but even with upsampling,
    grid resolution isn't infinite).
    Also for determining whether additional X-points are on the secondary
    separatrix after it is identified.
within_limiter_only:
    Only look for X-points that are within the limiting surface
"""
function x_points(
    g::GEQDSKFile;
    neighborhood::Int64=3,
    bpol_thresh_factor::Float64=20.0,
    up_factor::Int64=5,
    psin_tolerance::Float64=0.002,
    within_limiter_only::Bool=true,
)
    # Get out the basics
    r_eq = collect(g.r)
    z_eq = collect(g.z)
    psi_eq = g.psirz
    psia = g.simag
    psib = g.sibry
    psin_eq = (psi_eq .- psia) ./ (psib - psia)

    (br, bz, bpol, d) = fluxinfo(r_eq, z_eq, psi_eq)

    # Consider the wall because we don't care about field nulls outside of the VV
    wall_r = g.rlim
    wall_z = g.zlim
    wall_poly = [[wall_r[i], wall_z[i]] for i ∈ 1:length(wall_r)]
    wall_poly = append!(wall_poly, [[wall_r[1], wall_z[1]]])
    within_limiter =
        Matrix(
            reduce(
                hcat,
                [
                    [
                        PolygonOps.inpolygon([r_eq[i], z_eq[j]], wall_poly) for
                        j ∈ 1:length(z_eq)
                    ] for i ∈ 1:length(r_eq)
                ],
            )',
        ) .!= 0

    # Grid resolution is limited, so finding the exact minima may take extra work
    # The neighborhood should be big enough to cover the whole low spot around an
    # X-point and small enough that it doesn't catch multiple X-points. This is
    # tricky in the case of a snowflake. But snowflakes probably can't be formed
    # with reactor relevant coil sets, so forget them.
    local_minimum = find_local_minimum(bpol, neighborhood)

    # Identify a search region
    # due to limited grid resolution, not all X-points will result in the same
    # minimum value on the EFIT grid
    bpol_thresh = minimum(bpol) * bpol_thresh_factor
    # These guys are all probably the closest grid cells to X-points
    prime_real_estate =
        (bpol .< bpol_thresh) .& (d .< 0.0) .& (bpol .== local_minimum)
    if within_limiter_only
        prime_real_estate .&= within_limiter
    end
    search_centers = findall(prime_real_estate)

    # Log found X-points so they can be sorted by psi_N to identify whether they are
    # on the primary separatrix, secondary separatix, or even farther out. It is
    # unlikely that a third separatrix would have an X-point within the bore of the
    # machine, but it's not impossible.
    xrs = []
    xzs = []
    xpsins = []

    # For each search center, up sample the local grid to improve accuracy of the
    # final search result and record the X-point positions
    psi_int = Interpolations.linear_interpolation((r_eq, z_eq), psi_eq)
    for sc ∈ search_centers
        # Would this way have a significant speed advantage?
        # if PolygonOps.inpolygon([r_eq[sc[1]], z_eq[sc[2]]], wall_poly) != 1
        #     println(r_eq[sc[1]], ", ", z_eq[sc[2]], " is not within the limiter")
        #     continue
        # end
        imin = maximum([sc[1] - neighborhood, 1])
        imax = minimum([sc[1] + neighborhood, length(r_eq)])
        jmin = maximum([sc[2] - neighborhood, 1])
        jmax = minimum([sc[2] + neighborhood, length(z_eq)])
        sample = psin_eq[imin:imax, jmin:jmax]
        newi = (imax - imin) * up_factor + 1
        newj = (jmax - jmin) * up_factor + 1
        up_r =
            collect(LinRange(r_eq[imin], r_eq[imax], (imax - imin) * up_factor + 1))
        up_z =
            collect(LinRange(z_eq[jmin], z_eq[jmax], (jmax - jmin) * up_factor + 1))
        up_sample = psi_int.(up_r .+ (0.0 .* up_z)', up_z' .+ up_r .* 0.0)
        _, _, up_bpol, up_d = fluxinfo(up_r, up_z, up_sample)
        _, m = findmin(up_bpol)
        xr = up_r[m[1]]
        xz = up_z[m[2]]
        xpsin = (up_sample[m[1], m[2]] - psia) / (psib - psia)
        append!(xrs, xr)
        append!(xzs, xz)
        append!(xpsins, xpsin)

        # println("found x point at ", xr, ", ", xz, " with psiN = ", xpsin)
    end

    # Categorize the X-points; x-points on the primary or secondary separatrices
    # can be entered in special places
    if any(xpsins .> (1 + psin_tolerance))
        secondary_psin = minimum(xpsins[xpsins.>(1+psin_tolerance)])
    else
        secondary_psin = 10.0
    end
    problematic = xpsins .< (1 - psin_tolerance)
    primary = (1 - psin_tolerance) .<= xpsins .<= (1 + psin_tolerance)
    secondary = secondary_psin .<= xpsins .<= (secondary_psin + 2 * psin_tolerance)
    outer = xpsins .> (secondary_psin + 2 * psin_tolerance)
    xseps = primary .+ 2 .* secondary .+ 4 .* outer

    return xrs, xzs, xpsins, xseps
end
