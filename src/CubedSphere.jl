module CubedSphere

export sn, cn, conformal_cubed_sphere_mapping, conformal_cubed_sphere_inverse_mapping, cartesian_to_lat_lon
export compute_deviation_from_isotropy
export conformal_cubed_sphere_coordinates, optimized_non_uniform_conformal_cubed_sphere_coordinates

using Printf
using TaylorSeries

include("spherical_geometry.jl")
include("rancic_taylor_coefficients.jl")
include("conformal_cubed_sphere.jl")
include("generate_non_uniform_conformal_mapping_coordinates.jl")


end # module
