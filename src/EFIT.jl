module EFIT

include("io.jl")
export readg, writeg, GEQDSKFile

include("geqdsk_imas.jl")

include("utils.jl")
export triangularity, ellipticity, elongation, minor_radius, major_radius, aspect_ratio, elevation
export x_points, fluxinfo

const test_gfile = (@__DIR__)*"/../test/g000001.01000"
const test_gfile2 = (@__DIR__)*"/../test/g184833.03600"

end # module