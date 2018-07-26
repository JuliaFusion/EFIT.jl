__precompile__()

module EFIT

include("io.jl")
export readg, GEQDSKFile

include("utils.jl")
export triangularity, ellipticity, elongation, minor_radius, major_radius, aspect_ratio

end # module
