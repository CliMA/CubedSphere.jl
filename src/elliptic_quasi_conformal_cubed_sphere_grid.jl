using LinearAlgebra
using Statistics: mean, norm
using Random
using Rotations
using CairoMakie
using CubedSphere
using Dates
using Printf


function compute_deviation_from_orthogonality_at_point(X, Y, i, j)

    x1 = X[i, j]
    y1 = Y[i, j]
    x2 = X[i+1, j]
    y2 = Y[i+1, j]
    x3 = X[i+1, j+1]
    y3 = Y[i+1, j+1]
    x4 = X[i, j+1]
    y4 = Y[i, j+1]

    # Compute the vectors.
    v1 = [x2 - x1, y2 - y1]
    v2 = [x3 - x2, y3 - y2]
    v3 = [x4 - x3, y4 - y3]
    v4 = [x1 - x4, y1 - y4]

    # Compute the angles.
    angle1 = abs(acos(dot(normalize(v1), normalize(-v4))))
    angle2 = abs(acos(dot(normalize(v2), normalize(-v1))))
    angle3 = abs(acos(dot(normalize(v3), normalize(-v2))))
    angle4 = abs(acos(dot(normalize(v4), normalize(-v3))))

    # Compute the deviation from orthogonality.
    deviation_from_orthogonality_at_point = (abs(angle1 - π/2) + abs(angle2 - π/2) + abs(angle3 - π/2)
                                             + abs(angle4 - π/2))

    return deviation_from_orthogonality_at_point

end


function compute_deviation_from_orthogonality(X, Y)

    Nx, Ny = size(X)
    deviation_from_orthogonality = zeros(Nx, Ny)

    for i in 1:Nx-1
        for j in 1:Ny-1
            deviation_from_orthogonality[i, j] = compute_deviation_from_orthogonality_at_point(X, Y, i, j)
        end
    end

    return deviation_from_orthogonality

end


function compute_deviation_from_isotropy_at_point(X, Y, i, j)

    x1 = X[i, j]
    y1 = Y[i, j]
    x2 = X[i+1, j]
    y2 = Y[i+1, j]
    x3 = X[i+1, j+1]
    y3 = Y[i+1, j+1]
    x4 = X[i, j+1]
    y4 = Y[i, j+1]

    # Compute the vectors.
    v1 = [x2 - x1, y2 - y1]
    v2 = [x3 - x2, y3 - y2]
    v3 = [x4 - x3, y4 - y3]
    v4 = [x1 - x4, y1 - y4]

    # Compute the deviation from isotropy.
    deviation_from_isotropy_at_point = (
    abs(norm(v1) - norm(v2)) + abs(norm(v2) - norm(v3)) + abs(norm(v3) - norm(v4)) + abs(norm(v4) - norm(v1)))

    return deviation_from_isotropy_at_point

end


function compute_deviation_from_isotropy(X, Y)

    Nx, Ny = size(X)
    deviation_from_isotropy = zeros(Nx, Ny)

    for i in 1:Nx-1
        for j in 1:Ny-1
            deviation_from_isotropy[i, j] = compute_deviation_from_isotropy_at_point(X, Y, i, j)
        end
    end

    return deviation_from_isotropy

end


function compute_cell_area(X, Y, i, j)

    x1 = X[i,   j]
    y1 = Y[i,   j]
    x2 = X[i+1, j]
    y2 = Y[i+1, j]
    x3 = X[i+1, j+1]
    y3 = Y[i+1, j+1]
    x4 = X[i,   j+1]
    y4 = Y[i,   j+1]

    # Compute vectors.
    v1 = [x2 - x1, y2 - y1]
    v2 = [x3 - x2, y3 - y2]
    v3 = [x4 - x3, y4 - y3]
    v4 = [x1 - x4, y1 - y4]

    # Compute angles.
    angle1 = abs(acos(dot(normalize(v1), normalize(-v4))))
    angle3 = abs(acos(dot(normalize(v3), normalize(-v2))))

    # Compute area using triangles formed by the cell.
    s1 = 0.5 * norm(v1) * norm(v4) * sin(angle1)
    s3 = 0.5 * norm(v3) * norm(v2) * sin(angle3)

    cell_area = s1 + s3

    return cell_area

end


function compute_cell_areas(X, Y)

    Nx, Ny = size(X)
    cell_areas = zeros(Nx, Ny)

    for i in 1:Nx-1
        for j in 1:Ny-1
            cell_areas[i, j] = compute_cell_area(X, Y, i, j)
        end
    end

    return cell_areas

end


function conformal_cubed_sphere_coordinates(Nx, Ny)
    
    x = range(-1, 1, length = Nx)
    y = range(-1, 1, length = Ny)
    
    X = zeros(length(x), length(y))
    Y = zeros(length(x), length(y))
    Z = zeros(length(x), length(y))
    
    for (j, y′) in enumerate(y), (i, x′) in enumerate(x)
        X[i, j], Y[i, j], Z[i, j] = conformal_cubed_sphere_mapping(x′, y′)
    end
    
    return x, y, X, Y, Z
    
end


function derivative(y, x)

    N = length(x)
    dydx = zeros(N)
    
    for i in 1:N
        if i == 1
            dydx[i] = (y[i+1] - y[i])/(x[i+1] - x[i])
        elseif i == N
            dydx[i] = (y[i] - y[i-1])/(x[i] - x[i-1])
        else
            dydx[i] = (y[i+1] - y[i-1])/(x[i+1] - x[i-1])
        end
    end
    
end


function derivatives(u, x, y)

    Nx, Ny = size(u)

    dudx = zeros(Nx, Ny)
    dudy = zeros(Nx, Ny)
    
    for j in 1:Ny
        for i in 1:Nx
            if i == 1
                dudx[i, j] = (u[i+1, j] - u[i, j])/(x[i+1] - x[i])
            elseif i == Nx
                dudx[i, j] = (u[i, j] - u[i-1, j])/(x[i] - x[i-1])
            else
                dudx[i, j] = (u[i+1, j] - u[i-1, j])/(x[i+1] - x[i-1])
            end
        end
    end
    
    for j in 1:Ny
        for i in 1:Nx
            if j == 1
                dudy[i, j] = (u[i, j+1] - u[i, j])/(y[j+1] - y[j])
            elseif j == Ny
                dudy[i, j] = (u[i, j] - u[i, j-1])/(y[j] - y[j-1])
            else
                dudy[i, j] = (u[i, j+1] - u[i, j-1])/(y[j+1] - y[j-1])
            end
        end
    end
    
    return dudx, dudy
    
end


function specify_orthogonality_enforcing_control_function_components!(a_P, a_Q, b_P, b_Q, X_x, X_y, Y_x, Y_y, X_xx,
                                                                      X_xy, X_yy, Y_xx, Y_xy, Y_yy, S_x, S_y,
                                                                      boundary_identifiers, boundary_type,
                                                                      initialization, ξ, η, P, Q)
    
    ω_P = 0.05
    ω_Q = 0.05
    
    south = boundary_identifiers[1]
    east  = boundary_identifiers[2]
    north = boundary_identifiers[3]
    west  = boundary_identifiers[4]

    η₀ = (boundary_type == south) ? -1 : ((boundary_type == north) ?  1 : 0)
    ξ₀ = (boundary_type == east)  ?  1 : ((boundary_type == west)  ? -1 : 0)

    α = X_y^2 + Y_y^2
    β = X_x * X_y + Y_x * Y_y
    γ = X_x^2 + Y_x^2
    J = X_x * Y_y - Y_x * X_y
    
    R_1 = α * X_xx - 2β * X_xy + γ * X_yy
    R_2 = α * Y_xx - 2β * Y_xy + γ * Y_yy
    
    R = [R_1, R_2]
    
    if boundary_type == south || boundary_type == north
        X_y = - S_y * Y_x / sqrt(γ)
        Y_y =   S_y * X_x / sqrt(γ)
        a₁₁ = exp(-b_P * η₀) * X_x 
        a₁₂ = exp(-b_Q * η₀) * X_y
        a₂₁ = exp(-b_P * η₀) * Y_x
        a₂₂ = exp(-b_Q * η₀) * Y_y
    elseif boundary_type == east || boundary_type == west
        X_x = - S_x * Y_y / sqrt(α)
        Y_x =   S_x * X_y / sqrt(α)
        a₁₁ = exp(-b_P * ξ₀) * X_x 
        a₁₂ = exp(-b_Q * ξ₀) * X_y
        a₂₁ = exp(-b_P * ξ₀) * Y_x
        a₂₂ = exp(-b_Q * ξ₀) * Y_y
    end
    
    A = -J^2 * [a₁₁ a₁₂; a₂₁ a₂₂]
    PQ = A \ R

    if initialization
        P = PQ[1]
        Q = PQ[2]
    else
        P = P + ω_P * (PQ[1] - P)
        Q = Q + ω_Q * (PQ[2] - Q)
    end
    
    if boundary_type == south || boundary_type == north
        P = a_P * P * exp(-b_P * η)
        Q = a_Q * Q * exp(-b_Q * η)
    elseif boundary_type == east || boundary_type == west
        P = a_P * P * exp(-b_P * ξ)
        Q = a_Q * Q * exp(-b_Q * ξ)
    end

end 


function elliptic_quasi_conformal_cubed_sphere_coordinates(Nx, Ny, θ; return_model_diagnostics = false,
                                                           display_outcome = true)
    
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
    
    θ_x = θ[1]
    θ_y = θ[2]
    b_x = θ[3] # The decay scale, R_x = 1/b_x. So smaller b_x means larger decay scale R_x.
    b_y = θ[4] # The decay scale, R_y = 1/b_y. So smaller b_y means larger decay scale R_y.
    S_0 = θ[5]
    a_P = θ[6]
    a_Q = θ[7]
    b_P = θ[8]
    b_Q = θ[9]
    
    Δξ = 2/(Nx - 1)
    Δη = 2/(Ny - 1)
    
    a_x = θ_x * Δξ
    a_y = θ_y * Δη
    
    xArray = collect(x)
    yArray = collect(y)
    
    ΔxStretched = zeros(eltype(xArray), size(xArray))
    ΔyStretched = zeros(eltype(yArray), size(yArray))
    
    if Nx % 2 == 0
        Nx_mid = Nx ÷ 2
    else
        Nx_mid = (Nx + 1) ÷ 2 - 1
    end
    
    if Ny % 2 == 0
        Ny_mid = Ny ÷ 2
    else
        Ny_mid = (Ny + 1) ÷ 2 - 1
    end
    
    for i = 2:Nx_mid
        ΔxStretched[i] = a_x * exp(-b_x * (x[i] - x[2]))
        ΔxStretched[Nx - i + 1] = -ΔxStretched[i]
    end
    
    for j = 2:Ny_mid
        ΔyStretched[j] = a_y * exp(-b_y * (y[j] - y[2]))
        ΔyStretched[Ny - j + 1] = -ΔyStretched[j]
    end
    
    xStretched = x + ΔxStretched
    yStretched = y + ΔyStretched
    
    # Southern boundary    
    for (i, x′) in enumerate(xStretched)
        X[i, 1], Y[i, 1], Z[i, 1] = conformal_cubed_sphere_mapping(x′, y[1])
    end
    
    # Eastern boundary
    for (j, y′) in enumerate(yStretched)
        X[Nx, j], Y[Nx, j], Z[Nx, j] = conformal_cubed_sphere_mapping(x[Nx], y′)
    end
    
    # Northern boundary
    for (i, x′) in enumerate(xStretched)
        X[i, Ny], Y[i, Ny], Z[i, Ny] = conformal_cubed_sphere_mapping(x′, y[Ny])
    end
    
    # Western boundary
    for (j, y′) in enumerate(yStretched)
        X[1, j], Y[1, j], Z[1, j] = conformal_cubed_sphere_mapping(x[1], y′)
    end
    
    boundary_types = zeros(Nx, Ny)
    
    south = 1
    east  = 2
    north = 3
    west  = 4
    
    boundary_identifiers = [south, east, north, west]
    
    for j = 2:2
        for i = 2:Nx-1
            boundary_types[i, j] = south
        end
    end
    
    for j = 2:Ny-1
        for i = Nx-1:Nx-1
            boundary_types[i, j] = east
        end
    end
    
    for j = Ny-1:Ny-1
        for i = 2:Nx-1
            boundary_types[i, j] = north
        end
    end
    
    for j = 2:Ny-1
        for i = 2:2
            boundary_types[i, j] = west
        end
    end
    
    radius = 1
    S_x = S_0 * (sqrt(4π * radius^2/6))/2
    S_y = S_x
        
    P = zeros(Nx, Ny)
    Q = zeros(Nx, Ny)
    
    ω_x = 1.5
    ω_y = 1.5
    
    Residual_X = zeros(Nx, Ny)
    Residual_Y = zeros(Nx, Ny)
    
    nIterations = 10^6
    iteration_final = nIterations
    tolerance = 10^(-4)
    converged = false
    
    time_start = now()
    
    for iteration in 1:nIterations
    
        X_x, X_y = derivatives(X, x, y) 
        Y_x, Y_y = derivatives(Y, x, y)
        X_xx, X_xy = derivatives(X_x, x, y)
        X_xy, X_yy = derivatives(X_y, x, y)
        Y_xx, Y_xy = derivatives(Y_x, x, y)
        Y_xy, Y_yy = derivatives(Y_y, x, y)
        
        (iteration == 1) ? initialization = true : initialization = false
        
        for j = 2:Ny-1
            for i = 2:Nx-1
            
                if ((boundary_types[i, j] == south) || (boundary_types[i, j] == east) || (boundary_types[i, j] == north) 
                    || (boundary_types[i, j] == west))
                        specify_orthogonality_enforcing_control_function_components!(
                        a_P, a_Q, b_P, b_Q, X_x[i, j], X_y[i, j], Y_x[i, j], Y_y[i, j], X_xx[i, j], X_xy[i, j],
                        X_yy[i, j], Y_xx[i, j], Y_xy[i, j], Y_yy[i, j], S_x, S_y, boundary_identifiers,
                        boundary_types[i, j], initialization, x[i], y[j], P[i, j], Q[i, j])
                end
                
                α = X_y[i, j]^2 + Y_y[i, j]^2
                β = X_x[i, j] * X_y[i, j] + Y_x[i, j] * Y_y[i, j]
                γ = X_x[i, j]^2 + Y_x[i, j]^2
                J = X_x[i, j] * Y_y[i, j] - Y_x[i, j] * X_y[i, j]
                
                R_1 = α * X_xx[i, j] - 2β * X_xy[i, j] + γ * X_yy[i, j]
                R_2 = α * Y_xx[i, j] - 2β * Y_xy[i, j] + γ * Y_yy[i, j]
                
                RHS_X = -J^2 * (P[i, j] * X_x[i, j] + Q[i, j] * X_y[i, j])
                RHS_Y = -J^2 * (P[i, j] * Y_x[i, j] + Q[i, j] * Y_y[i, j])
                
                Residual_X[i, j] = RHS_X - R_1
                Residual_Y[i, j] = RHS_Y - R_2
                
                X[i, j] = X[i, j] - ω_x * Residual_X[i, j]/(2(α/Δξ^2 + γ/Δη^2))
                Y[i, j] = Y[i, j] - ω_y * Residual_Y[i, j]/(2(α/Δξ^2 + γ/Δη^2))
                
            end
        end
        
        Residual_maximum = max(maximum(abs, Residual_X), maximum(abs, Residual_Y))
        
        if Residual_maximum < tolerance
            converged = true
            iteration_final = iteration
            break
        end
    
    end
    
    for j in 1:Ny
        for i in 1:Nx
            Z[i, j] = sqrt(1 - X[i, j]^2 - Y[i, j]^2)
        end
    end
    
    time_end = now()
    time_elapsed = time_end - time_start
    
    if display_outcome
        if converged
            @printf("The SOR has converged after %d iterations.\n", iteration_final)
            @printf("The number of iterations required for convergence is %d.\n", iteration_final)
            @printf("Elapsed CPU time is %s.\n", time_elapsed)
        else
            @printf("The SOR has not converged after %d iterations.\n", iteration_final)
        end
        @printf("Average time per iteration is %.3g.\n", Dates.value(time_elapsed)/(1000 * iteration_final))
    end

    if return_model_diagnostics

        deviation_from_orthogonality = compute_deviation_from_orthogonality(X, Y)
        deviation_from_isotropy = compute_deviation_from_isotropy(X, Y)
        cell_areas = compute_cell_areas(X, Y)
        minimum_to_maximum_cell_area_ratio = minimum(cell_areas[1:Nx-1,1:Ny-1])/maximum(cell_areas[1:Nx-1,1:Ny-1])

        model_diagnostics = vcat(deviation_from_orthogonality[:], deviation_from_isotropy[:],
                                 minimum_to_maximum_cell_area_ratio)

        return x, y, X, Y, Z, model_diagnostics

    else
    
        return x, y, X, Y, Z

    end
    
end


function visualize_conformal_cubed_sphere_2D_3D(X, Y, Z, filename)

    axis_kwargs_2D = (xlabelsize = 22.5, ylabelsize = 22.5, xticklabelsize = 17.5, yticklabelsize = 17.5, 
                      xticklabelpad = 10, yticklabelpad = 10, titlesize = 27.5, titlegap = 15, titlefont = :bold, 
                      xlabel = "x", ylabel = "y")
    axis_kwargs_3D = (xlabelsize = 22.5, ylabelsize = 22.5, zlabelsize = 22.5, xticklabelsize = 17.5, 
                      yticklabelsize = 17.5, zticklabelsize = 17.5, xticklabelpad = 10, yticklabelpad = 10, 
                      zticklabelpad = 10, titlesize = 27.5, titlegap = 15, titlefont = :bold, xlabel = "x", 
                      ylabel = "y", zlabel = "z")
    hide_decorations = false

    colors = [:orange, :red, :deepskyblue, :purple, :green, :blue]
    alphas = [1, 1, 0.1125, 0.1125, 1, 0.1125]

    fig = Figure(resolution = (1600, 750))

    ax2D = Axis(fig[1, 1]; aspect = 1, title = "Conformal Cubed Sphere Projection in 2D", axis_kwargs_2D...)
    ax3D = Axis3(fig[1, 2]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)),
                 title = "Conformal Cubed Sphere in 3D", axis_kwargs_3D...)
    
    for ax in [ax2D, ax3D]
        hide_decorations && hidedecorations!(ax)
        wireframe!(ax, X, Y, Z, color = colors[1], alpha = alphas[1])
    end
    
    rotations = (RotX(π/2), RotX(-π/2), RotY(π/2), RotY(-π/2), RotX(π))
    
    for (i, R) in enumerate(rotations)
    
        X′ = similar(X)
        Y′ = similar(Y)
        Z′ = similar(Z)
    
        for I in CartesianIndices(X)
            X′[I], Y′[I], Z′[I] = R * [X[I], Y[I], Z[I]]
        end
        
        wireframe!(ax2D, X′, Y′, Z′, color = colors[i+1])
        wireframe!(ax3D, X′, Y′, Z′, color = colors[i+1], alpha = alphas[i+1])
        
    end
    
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 200)
    
    save(filename, fig)
    
end


function specify_parameters()

    θ_x = 0.75
    θ_y = 0.75
    b_x = 3
    b_y = 3
    S_0 = 1
    a_P = 10^(-4)
    a_Q = 10^(-4)
    b_P = 0.5
    b_Q = 0.5

    θ = [θ_x, θ_y, b_x, b_y, S_0, a_P, a_Q, b_P, b_Q]

    return θ

end


function specify_parameter_limits()

    θ_x_limits = [0.5, 1.25]
    θ_y_limits = [0.5, 1.25]
    b_x_limits = [2, 4]
    b_y_limits = [2, 4]
    S_0_limits = [0.25, 2.5]
    a_P_limits = [10^(-4), 10^(-1)]
    a_Q_limits = [10^(-4), 10^(-1)]
    b_P_limits = [0.1, 1]
    b_Q_limits = [0.1, 1]

    θ_limits = [θ_x_limits, θ_y_limits, b_x_limits, b_y_limits, S_0_limits, a_P_limits, a_Q_limits, b_P_limits,
                b_Q_limits]

    return θ_limits

end


function enforce_symmetry_in_parameters!(θ)

    θ[2] = θ[1]
    θ[4] = θ[3]
    θ[7] = θ[6]
    θ[9] = θ[8]

end


function specify_random_parameters(nEnsemble; enforce_symmetry = false)

    θ = specify_parameters()
    θ_limits = specify_parameter_limits()

    θᵣ = [[θ_limits[j][1] + (θ_limits[j][2] - θ_limits[j][1]) * rand() for j in 1:lastindex(θ)] for i in 1:nEnsemble]

    if enforce_symmetry
        for i in 1:nEnsemble
            enforce_symmetry_in_parameters!(θᵣ[i])
        end
    end

    return θᵣ

end


function visualize_conformal_and_elliptic_quasi_conformal_cubed_sphere_2D_3D()

    Nx, Ny = 32 + 1, 32 + 1
    θ = specify_parameters()

    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
    filename = "conformal_cubed_sphere_2D_3D.png"
    visualize_conformal_cubed_sphere_2D_3D(X, Y, Z, filename)

    cell_areas = compute_cell_areas(X, Y)
    minimum_to_maximum_cell_area_ratio = minimum(cell_areas[1:Nx-1,1:Ny-1])/maximum(cell_areas[1:Nx-1,1:Ny-1])

    x, y, X, Y, Z = elliptic_quasi_conformal_cubed_sphere_coordinates(Nx, Ny, θ)
    filename = "elliptic_quasi_conformal_cubed_sphere_2D_3D.png"
    visualize_conformal_cubed_sphere_2D_3D(X, Y, Z, filename)

    cell_areas = compute_cell_areas(X, Y)
    minimum_to_maximum_cell_area_ratio = minimum(cell_areas[1:Nx-1,1:Ny-1])/maximum(cell_areas[1:Nx-1,1:Ny-1])

end


visualize_conformal_and_elliptic_quasi_conformal_cubed_sphere_2D_3D()


function specify_weights_for_model_diagnostics()

    weights = [10, 1, 10]

    return weights

end


function compute_weighted_model_diagnostics(Nx, Ny, model_diagnostics)

    deviation_from_orthogonality = model_diagnostics[1:Nx*Ny]
    deviation_from_isotropy = model_diagnostics[Nx*Ny+1:2*Nx*Ny]
    minimum_to_maximum_cell_area_ratio = model_diagnostics[2*Nx*Ny+1]

    weights = specify_weights_for_model_diagnostics()

    weighted_model_diagnostics = vcat(weights[1] * deviation_from_orthogonality[:],
                                      weights[2] * deviation_from_isotropy[:],
                                      weights[3] * minimum_to_maximum_cell_area_ratio)

    return weighted_model_diagnostics

end


function forward_map(Nx, Ny, θ)

    θ_limits = specify_parameter_limits()

    for i in 1:lastindex(θ)
        θ[i] = clamp(θ[i], θ_limits[i][1], θ_limits[i][2])
    end

    if Nx == Ny
        enforce_symmetry_in_parameters!(θ)
    end

    x, y, X, Y, Z, model_diagnostics = (
    elliptic_quasi_conformal_cubed_sphere_coordinates(Nx, Ny, θ; return_model_diagnostics = true,
                                                      display_outcome = false))

    weighted_model_diagnostics = compute_weighted_model_diagnostics(Nx, Ny, model_diagnostics)

	return weighted_model_diagnostics

end


function specify_ideal_weighted_model_diagnostics(Nx, Ny)

    deviation_from_orthogonality = zeros(Nx, Ny)
    deviation_from_isotropy = zeros(Nx, Ny)
    minimum_to_maximum_cell_area_ratio = 1

    weights = specify_weights_for_model_diagnostics()

    ideal_weighted_model_diagnostics = vcat(weights[1] * deviation_from_orthogonality[:],
                                            weights[2] * deviation_from_isotropy[:],
                                            weights[3] * minimum_to_maximum_cell_area_ratio)

    return ideal_weighted_model_diagnostics

end


function optimize!(Nx, Ny, θ; nIterations = 10, Δt = 1)

    ideal_data = specify_ideal_weighted_model_diagnostics(Nx, Ny)
    model_data = forward_map(Nx, Ny, mean(θ))

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

        #=
		Evaluating the forward map for all ensemble members. This is the most expensive step because it need to  run the
        model nEnsemble times. For the moment our model is simple, but imagine doing this with a full climate model!
		Luckily this step is embarassingly parallelizeable.
        =#
        Threads.@threads for n in 1:nEnsemble
			G[n] .= forward_map(Nx, Ny, θ[n]) # Error handling needs to go here.
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


function optimize_and_visualize_elliptic_quasi_conformal_cubed_sphere_2D_3D()

    Nx, Ny = 32 + 1, 32 + 1

    Nx == Ny ? enforce_symmetry = true : enforce_symmetry = false

    nEnsemble = 36 # Choose nEnsemble to be at least 4 times the number of parameters.

    begin

        Random.seed!(123)
        θᵣ = specify_random_parameters(nEnsemble; enforce_symmetry = enforce_symmetry)
        θᵢ = deepcopy(θᵣ)

        θ_series = optimize!(Nx, Ny, θᵣ; nIterations = 10)

    end

    @printf("\nThe unoptimized parameters are:\n\nθ_x = %.2f\nθ_y = %.2f\nb_x = %.2f\nb_y = %.2f\nS_0 = %.2f\na_P = %.2f\na_Q = %.2f\nb_P = %.2f\nb_Q = %.2f\n",
            mean(θᵢ)[1], mean(θᵢ)[2], mean(θᵢ)[3], mean(θᵢ)[4], mean(θᵢ)[5], mean(θᵢ)[6], mean(θᵢ)[7], mean(θᵢ)[8], mean(θᵢ)[9])

    @printf("\nThe optimized parameters are:\n\nθ_x = %.2f\nθ_y = %.2f\nb_x = %.2f\nb_y = %.2f\nS_0 = %.2f\na_P = %.2f\na_Q = %.2f\nb_P = %.2f\nb_Q = %.2f\n\n",
            mean(θᵣ)[1], mean(θᵣ)[2], mean(θᵣ)[3], mean(θᵣ)[4], mean(θᵣ)[5], mean(θᵣ)[6], mean(θᵣ)[7], mean(θᵣ)[8], mean(θᵣ)[9])

    x, y, X, Y, Z = elliptic_quasi_conformal_cubed_sphere_coordinates(Nx, Ny, mean(θᵣ))

    filename = "optimized_elliptic_quasi_conformal_cubed_sphere_2D_3D.png"
    visualize_conformal_cubed_sphere_2D_3D(X, Y, Z, filename)

end


optimize_and_visualize_elliptic_quasi_conformal_cubed_sphere_2D_3D()