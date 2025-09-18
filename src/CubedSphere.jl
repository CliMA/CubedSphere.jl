module CubedSphere

export sn, cn, conformal_cubed_sphere_mapping, conformal_cubed_sphere_inverse_mapping, cartesian_to_lat_lon
export lat_lon_to_cartesian, cartesian_to_lat_lon, turning_angle, spherical_distance, spherical_area_triangle,
    spherical_area_quadrilateral, spherical_quadrilateral_vertices, compute_deviation_from_isotropy, compute_cell_areas
export conformal_cubed_sphere_coordinates, optimized_non_uniform_conformal_cubed_sphere_coordinates

using Printf
using TaylorSeries

include("rancic_taylor_coefficients.jl")
include("conformal_cubed_sphere.jl")
include("spherical_geometry.jl")
include("generate_non_uniform_conformal_mapping_coordinates.jl")

end # module
