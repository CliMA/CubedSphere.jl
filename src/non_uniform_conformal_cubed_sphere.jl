using LinearAlgebra
using Statistics
using Random


"""
    spherical_distance(a₁::AbstractVector, a₂::AbstractVector)

Compute the arc length (distance) between two points on the unit sphere, given their Cartesian coordinates `a₁` and
`a₂`. The points `a₁` and `a₂` should be unit vectors (magnitude 1), representing locations on the sphere.
"""
function spherical_distance(a₁::AbstractVector, a₂::AbstractVector)
    (sum(a₁.^2) ≈ 1 && sum(a₂.^2) ≈ 1) || error("a₁ and a₂ must be unit vectors")

    # Compute the dot product and calculate the arccosine to find the angle.
    cosθ = dot(a₁, a₂)

    # Ensure the result is within the domain of acos due to potential floating-point errors.
    cosθ = clamp(cosθ, -1.0, 1.0)

    # Return the arc length, which is the angle between the two points.
    return acos(cosθ)
end


"""
    spherical_area_triangle(a::AbstractVector, b::AbstractVector, c::AbstractVector)

Return the area of a spherical triangle on the unit sphere with vertices given by the 3-vectors `a`, `b`, and `c` whose 
origin is the the center of the sphere. The formula was first given by Eriksson (1990).

If we denote with ``A``, ``B``, and ``C`` the inner angles of the spherical triangle and with ``a``, ``b``, and ``c``
the side of the triangle then, it has been known since Euler and Lagrange that
``\\tan(E/2) = P / (1 + \\cos a + \\cos b + \\cos c)``, where ``E = A + B + C - π`` is the triangle's excess and
``P = (1 - \\cos²a - \\cos²b - \\cos²c + 2 \\cos a \\cos b \\cos c)^{1/2}``.

On the unit sphere, ``E`` is precisely the area of the spherical triangle. Erikkson (1990) showed that ``P`` above is
the same as the volume defined by the vectors `a`, `b`, and `c`, that is ``P = |𝐚 \\cdot (𝐛 \\times 𝐜)|``.

References
==========
* Eriksson, F. (1990) On the measure of solid angles, Mathematics Magazine, 63 (3), 184-187, doi:10.1080/0025570X.1990.11977515
"""
function spherical_area_triangle(a₁::AbstractVector, a₂::AbstractVector, a₃::AbstractVector)
    (sum(a₁.^2) ≈ 1 && sum(a₂.^2) ≈ 1 && sum(a₃.^2) ≈ 1) || error("a₁, a₂, a₃ must be unit vectors")

    tan½E = abs(dot(a₁, cross(a₂, a₃)))
    tan½E /= 1 + dot(a₁, a₂) + dot(a₂, a₃) + dot(a₁, a₃)

    return 2atan(tan½E)
end


"""
    spherical_area_quadrilateral(a₁, a₂, a₃, a₄)

Return the area of a spherical quadrilateral on the unit sphere whose points are given by 3-vectors, `a`, `b`, `c`, and
`d`. The area of the quadrilateral is given as the sum of the ares of the two non-overlapping triangles. To avoid having
to pick the triangles appropriately ensuring they are not overlapping, we compute the area of the quadrilateral as the
half the sum of the areas of all four potential triangles formed by `a₁`, `a₂`, `a₃`, and `a₄`.
"""
spherical_area_quadrilateral(a::AbstractVector, b::AbstractVector, c::AbstractVector, d::AbstractVector) =
    1/2 * (spherical_area_triangle(a, b, c) + spherical_area_triangle(a, b, d) +
           spherical_area_triangle(a, c, d) + spherical_area_triangle(b, c, d))


function spherical_quadrilateral_vertices(X, Y, Z, i, j)
    x₁ = X[i, j]
    y₁ = Y[i, j]
    z₁ = Z[i, j]
    a₁ = [x₁, y₁, z₁]
    x₂ = X[i+1, j]
    y₂ = Y[i+1, j]
    z₂ = Z[i+1, j]
    a₂ = [x₂, y₂, z₂]
    x₃ = X[i+1, j+1]
    y₃ = Y[i+1, j+1]
    z₃ = Z[i+1, j+1]
    a₃ = [x₃, y₃, z₃]
    x₄ = X[i, j+1]
    y₄ = Y[i, j+1]
    z₄ = Z[i, j+1]
    a₄ = [x₄, y₄, z₄]

    return a₁, a₂, a₃, a₄
end


function compute_deviation_from_isotropy(X, Y, Z)
    Nx, Ny = size(X)
    deviation_from_isotropy = zeros(Nx-1, Ny-1)

    for j in 1:Ny-1, i in 1:Nx-1
        a₁, a₂, a₃, a₄ = spherical_quadrilateral_vertices(X, Y, Z, i, j)

        # Compute the arc lengths (distances) between the points a₁ and a₂, a₂ and a₃, a₃ and a₄, and a₄ and a₁ on the
        # unit sphere.
        d₁ = spherical_distance(a₁, a₂)
        d₂ = spherical_distance(a₂, a₃)
        d₃ = spherical_distance(a₃, a₄)
        d₄ = spherical_distance(a₄, a₁)

        # Compute the deviation from isotropy.
        deviation_from_isotropy[i, j] = abs(d₁ - d₂) + abs(d₂ - d₃) + abs(d₃ - d₄) + abs(d₄ - d₁)
    end

    return norm(deviation_from_isotropy)
end


function compute_cell_areas(X, Y, Z)
    Nx, Ny = size(X)
    cell_areas = zeros(Nx-1, Ny-1)

    for j in 1:Ny-1, i in 1:Nx-1
        a₁, a₂, a₃, a₄ = spherical_quadrilateral_vertices(X, Y, Z, i, j)
        cell_areas[i, j] = spherical_area_quadrilateral(a₁, a₂, a₃, a₄)
    end

    return cell_areas
end


function geometric_spacing(N, ratio_raised_to_N_minus_one)
    ratio = ratio_raised_to_N_minus_one^(1/(N - 1))
    x_faces = zeros(N)

    if isodd(N)
        M = round(Int, (N + 1)/2)
    
        Δx = 1 * (ratio - 1) / (ratio^(M - 1) - 1)

        x_faces[M] = 0
        
        k = 0
        
        for i in M+1:N
            x_faces[i] = x_faces[i-1] + Δx * ratio^k
            x_faces[N+1-i] = -x_faces[i]
            k += 1
        end
        
        x_faces[1] = -1
        x_faces[N] = 1
    else
        M = Int(N/2)
    
        Δx = 1/((ratio^M - 1)/(ratio - 1) - 0.5)
        
        x_faces[M] = -0.5Δx
        x_faces[M+1] = 0.5Δx
        
        k = 1
        
        for i in M+2:N
            x_faces[i] = x_faces[i-1] + Δx * ratio^k
            x_faces[N+1-i] = -x_faces[i]
            k += 1
        end
        
        x_faces[1] = -1
        x_faces[N] = 1
    end 
    
    return x_faces
end


function exponential_spacing(N, k₀ByN)
    k₀ = k₀ByN * N
    x_faces = zeros(N)
    
    if isodd(N)
        M = round(Int, (N + 1)/2)
        
        A = [exp(1/k₀) 1
             exp(M/k₀) 1]

        b = [0, 1]
        
        coefficients = A \ b
        
        x_faces[M:N] = coefficients[1] * exp.((1:M)/k₀) .+ coefficients[2]
        
        for i in 1:M-1
            x_faces[i] = -x_faces[N+1-i]
        end
        
        x_faces[1] = -1
        x_faces[M] = 0
        x_faces[N] = 1
    else
        M = Int(N/2)
        
        A = [exp(1.5/k₀)   1
             exp((M+1)/k₀) 1]
    
        b = [0, 1]
        
        coefficients = A \ b
        
        x_faces[M+1:N] = coefficients[1] * exp.((2:M+1)/k₀) .+ coefficients[2]
        
        for i in 1:M
            x_faces[i] = -x_faces[N+1-i]
        end
        
        x_faces[1] = -1
        x_faces[N] = 1
    end

    return x_faces
end


function conformal_cubed_sphere_coordinates(Nx, Ny; non_uniform_spacing = false, spacing_type = "geometric",
                                            ratio_raised_to_Nx_minus_one = 10.5, k₀ByNx = 0.45)
    x = range(-1, 1, length = Nx)
    y = range(-1, 1, length = Ny)

    if non_uniform_spacing
        if spacing_type == "geometric"
            # For Nx = Ny = 32 + 1, setting ratio = 1.0775 increases the minimum cell width by a factor of 1.92.
            # For Nx = Ny = 1024 + 1, setting ratio = 1.0042 increases the minimum cell width by a factor of 3.25.
            x = geometric_spacing(Nx, ratio_raised_to_Nx_minus_one)
            y = geometric_spacing(Ny, ratio_raised_to_Nx_minus_one)
        elseif spacing_type == "exponential"
            # For Nx = Ny = 32 + 1, setting k₀ByNx = 15 increases the minimum cell width by a factor of 1.84.
            # For Nx = Ny = 1024 + 1, setting k₀ByNx = 10 increases the minimum cell width by a factor of 2.58.
            x = exponential_spacing(Nx, k₀ByNx)
            y = exponential_spacing(Ny, k₀ByNx)
        end
    end
    
    X = zeros(length(x), length(y))
    Y = zeros(length(x), length(y))
    Z = zeros(length(x), length(y))
    
    for (j, y′) in enumerate(y), (i, x′) in enumerate(x)
        X[i, j], Y[i, j], Z[i, j] = conformal_cubed_sphere_mapping(x′, y′)
    end
    
    return x, y, X, Y, Z
end


function specify_parameters(spacing_type)
    θ = 0.0
    if spacing_type == "geometric"
        ratio_raised_to_N = 1.0775
        θ = ratio_raised_to_N
    elseif spacing_type == "exponential"
        k₀ByN = 15
        θ = k₀ByN
    end
    return [θ]
end


function specify_parameter_limits(spacing_type)
    θ_limits = zeros(2)
    if spacing_type == "geometric"
        ratio_raised_to_N_limits = [5, 15]
        θ_limits[1] = ratio_raised_to_N_limits[1]
        θ_limits[2] = ratio_raised_to_N_limits[2]
    elseif spacing_type == "exponential"
        k₀ByN_limits = [0.4, 0.5]
        θ_limits[1] = k₀ByN_limits[1]
        θ_limits[2] = k₀ByN_limits[2]
    end
    return [θ_limits]
end


function specify_random_parameters(nEnsemble, spacing_type)
    θ = specify_parameters(spacing_type)
    θ_limits = specify_parameter_limits(spacing_type)

    θᵣ = [[θ_limits[j][1] + (θ_limits[j][2] - θ_limits[j][1]) * rand() for j in 1:lastindex(θ)] for i in 1:nEnsemble]

    return θᵣ
end


function specify_weights_for_model_diagnostics()
    weights = [10, 1]
    return weights
end


function compute_model_diagnostics(X, Y, Z, minimum_reference_cell_area)
    cell_areas = compute_cell_areas(X, Y, Z)
    normalized_minimum_cell_width = sqrt(minimum(cell_areas)/minimum_reference_cell_area)
    
    deviation_from_isotropy = compute_deviation_from_isotropy(X, Y, Z)

    model_diagnostics = vcat(normalized_minimum_cell_width, deviation_from_isotropy)
    
    return model_diagnostics
end


function compute_weighted_model_diagnostics(model_diagnostics)
    normalized_minimum_cell_width = model_diagnostics[1]
    deviation_from_isotropy = model_diagnostics[2]

    weights = specify_weights_for_model_diagnostics()

    weighted_model_diagnostics = vcat(weights[1] * normalized_minimum_cell_width, weights[2] * deviation_from_isotropy)

    return weighted_model_diagnostics
end


function forward_map(Nx, Ny, spacing_type, θ)
    θ_limits = specify_parameter_limits(spacing_type)

    for i in 1:lastindex(θ)
        θ[i] = clamp(θ[i], θ_limits[i][1], θ_limits[i][2])
    end
    
    x_reference, y_reference, X_reference, Y_reference, Z_reference = conformal_cubed_sphere_coordinates(Nx, Ny)
    cell_areas = compute_cell_areas(X_reference, Y_reference, Z_reference)
    minimum_reference_cell_area = minimum(cell_areas)
    
    x, y, X, Y, Z = (
    conformal_cubed_sphere_coordinates(Nx, Ny; non_uniform_spacing = true, spacing_type = spacing_type,
                                       ratio_raised_to_Nx_minus_one = θ[1], k₀ByNx = θ[1]))
    
    model_diagnostics = compute_model_diagnostics(X, Y, Z, minimum_reference_cell_area)

    weighted_model_diagnostics = compute_weighted_model_diagnostics(model_diagnostics)

	return weighted_model_diagnostics
end


function specify_ideal_weighted_model_diagnostics()
    normalized_minimum_cell_width = 4

    deviation_from_isotropy = 0.0

    weights = specify_weights_for_model_diagnostics()

    ideal_weighted_model_diagnostics = vcat(weights[1] * normalized_minimum_cell_width,
                                            weights[2] * deviation_from_isotropy)

    return ideal_weighted_model_diagnostics
end


function optimize!(Nx, Ny, spacing_type, θ; nIterations = 10, Δt = 1)
    ideal_data = specify_ideal_weighted_model_diagnostics()
    model_data = forward_map(Nx, Ny, spacing_type, mean(θ))

    nData = length(ideal_data)
    nEnsemble = length(θ)

    θ_series = [copy(θ)]

    error = norm(model_data - ideal_data)

    @printf("\n")
    @info "Iteration 0 with error $error"

	G = [copy(model_data) for i in 1:nEnsemble]

	# EKI iteration is equivalent to a time step of the above equation.
    @inbounds for i in 1:nIterations
        θ̄ = mean(θ)

		# Evaluating the forward map for all ensemble members. This is the most expensive step because it needs to run
        # the model nEnsemble times. For the moment our model is simple, but imagine doing this with a full climate
        # model! Luckily this step is embarassingly parallelizeable.
        
        Threads.@threads for n in 1:nEnsemble
			G[n] .= forward_map(Nx, Ny, spacing_type, θ[n]) # Error handling needs to go here.
		end

		# The ensemble mean output of the models
		G̅ = mean(G)

        # Calculating the covariances to be used in the update steps
        Cᵘᵖ = (θ[1] - θ̄) * (G[1] - G̅)'
        Cᵖᵖ = (G[1] - G̅) * (G[1] - G̅)'

        for j = 2:nEnsemble
            Cᵘᵖ += (θ[j] - θ̄) * (G[j] - G̅)'
            Cᵖᵖ += (G[j] - G̅) * (G[j] - G̅)'
        end

        Cᵘᵖ *= 1 / (nEnsemble - 1)
        Cᵖᵖ *= 1 / (nEnsemble - 1)

        # Ensemblize the data (adding the random noise η).
        y = [ideal_data + Δt * randn(nData) for i in 1:nEnsemble]

		# The residual from our observations
        r = y - G

        # Update the parameters using implicit pseudo-time-stepping, which involves solving a linear system.
        Cᵖᵖ_factorized = cholesky(Symmetric(Cᵖᵖ + 1 / Δt * LinearAlgebra.I))

        for j in 1:nEnsemble
            θ[j] .+= Cᵘᵖ * (Cᵖᵖ_factorized \ r[j])
        end

        error = norm(mean(r))
        @info "Iteration $i with error $error"
        push!(θ_series, copy(θ))
    end

    return θ_series
end


function optimized_non_uniform_conformal_cubed_sphere_coordinates(Nx, Ny, spacing_type)
    nEnsemble = 40 # Choose nEnsemble to be at least 4 times the number of parameters.
    
    @info "Optimize non-uniform conformal cubed sphere for Nx = $Nx and Ny = $Ny"
    
    begin
        Random.seed!(123)
        θᵣ = specify_random_parameters(nEnsemble, spacing_type)
        θᵢ = deepcopy(θᵣ)

        θ_series = optimize!(Nx, Ny, spacing_type, θᵣ; nIterations = 10)
    end
    
    if spacing_type == "geometric"
        θ_name = "ratio_raised_to_Nx_minus_one"
    elseif spacing_type == "exponential"
        θ_name = "k₀ByNx"
    end

    println("\nThe unoptimized parameters are: $θ_name = $(round(mean(θᵢ)[1], digits=2))\n")
    println("\nThe optimized parameters are: $θ_name = $(round(mean(θᵣ)[1], digits=2))\n")
    
    x, y, X, Y, Z = (
    conformal_cubed_sphere_coordinates(Nx, Ny; non_uniform_spacing = true, spacing_type = spacing_type,
                                       ratio_raised_to_Nx_minus_one = mean(θᵣ)[1], k₀ByNx = mean(θᵣ)[1]))
    
    return x, y, X, Y, Z
end
