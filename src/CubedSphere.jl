module CubedSphere

export sn, cn, conformal_cubed_sphere_mapping, conformal_cubed_sphere_inverse_mapping, cartesian_to_lat_lon
export compute_cell_areas, conformal_cubed_sphere_coordinates, optimized_non_uniform_conformal_cubed_sphere_coordinates

using Printf
using TaylorSeries

include("rancic_taylor_coefficients.jl")
include("conformal_cubed_sphere.jl")
include("cartesian_to_lat_lon.jl")
include("non_uniform_conformal_cubed_sphere.jl")

end # module
