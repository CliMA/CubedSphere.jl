module CubedSphere

export conformal_cubed_sphere_mapping, conformal_cubed_sphere_inverse_mapping, cartesian_to_lat_lon

using TaylorSeries

include("rancic_taylor_coefficients.jl")
include("conformal_cubed_sphere.jl")
include("cartesian_to_lat_lon.jl")

end # module
