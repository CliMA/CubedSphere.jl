module CubedSphere

export sn, cn

using Printf
using Requires
using TaylorSeries

include("complex_jacobi_elliptic.jl")
include("rancic_taylor_coefficients.jl")
include("conformal_cubed_sphere.jl")

function __init__()
    @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" begin
        include("visualize_conformal_mapping.jl")
    end
end

end # module
