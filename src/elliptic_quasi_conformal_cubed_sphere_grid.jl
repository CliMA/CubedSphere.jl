using Rotations
using CairoMakie
using CubedSphere
using Dates
using Printf


function conformal_cubed_sphere_coordinates(Nx, Ny, return_computational_coordinates = false)
    
    x = range(-1, 1, length = Nx)
    y = range(-1, 1, length = Ny)
    
    X = zeros(length(x), length(y))
    Y = zeros(length(x), length(y))
    Z = zeros(length(x), length(y))
    
    for (j, y′) in enumerate(y), (i, x′) in enumerate(x)
        X[i, j], Y[i, j], Z[i, j] = conformal_cubed_sphere_mapping(x′, y′)
    end
    
    if return_computational_coordinates
        return x, y, X, Y, Z
    else
        return X, Y, Z
    end
    
end


function derivative(y, x)

    N = length(x)
    dydx = zeros(N)
    
    for i = 1:N
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
    
    for j = 1:Ny
        for i = 1:Nx
            if i == 1
                dudx[i, j] = (u[i+1, j] - u[i, j])/(x[i+1] - x[i])
            elseif i == Nx
                dudx[i, j] = (u[i, j] - u[i-1, j])/(x[i] - x[i-1])
            else
                dudx[i, j] = (u[i+1, j] - u[i-1, j])/(x[i+1] - x[i-1])
            end
        end
    end
    
    for j = 1:Ny
        for i = 1:Nx
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


function specify_orthogonality_enforcing_control_function_components(X_x, X_y, Y_x, Y_y, X_xx, X_xy, X_yy, Y_xx, Y_xy, 
                                                                     Y_yy, S_x, S_y, boundary_identifiers, 
                                                                     boundary_type, initialization, ξ, η, P, Q)
    
    ω_P = 0.05
    ω_Q = 0.05
    
    a_P = 10^(-4)
    a_Q = 10^(-4)
    
    b_P = 0.5
    b_Q = 0.5
    
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
    
    return P, Q

end 


function elliptic_quasi_conformal_cubed_sphere_coordinates(Nx, Ny, return_computational_coordinates = false)
    
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny, true)
    
    θ_x = 0.75
    θ_y = 0.75
    
    Δξ = 2/(Nx - 1)
    Δη = 2/(Ny - 1)
    
    a_x = θ_x * Δξ
    a_y = θ_y * Δη
    
    b_x = 3 # The decay scale, R_x = 1/b_x. So smaller b_x means larger decay scale R_x.
    b_y = 3 # The decay scale, R_y = 1/b_y. So smaller b_y means larger decay scale R_y.
    
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
    S_x = (sqrt(4π * radius^2/6))/2
    S_y = S_x
        
    P = zeros(Nx, Ny)
    Q = zeros(Nx, Ny)
    
    ω_x = 1.5
    ω_y = 1.5
    
    Residual_X = zeros(Nx, Ny)
    Residual_Y = zeros(Nx, Ny)
    
    n_iterations = 10^6
    iteration_final = n_iterations
    tolerance = 10^(-4)
    converged = false
    
    time_start = now()
    
    for iteration = 1:n_iterations 
    
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
                        P[i, j], Q[i, j] = (
                        specify_orthogonality_enforcing_control_function_components(
                        X_x[i, j], X_y[i, j], Y_x[i, j], Y_y[i, j], X_xx[i, j], X_xy[i, j], X_yy[i, j], Y_xx[i, j], 
                        Y_xy[i, j], Y_yy[i, j], S_x, S_y, boundary_identifiers, boundary_types[i, j], initialization, 
                        x[i], y[j], P[i, j], Q[i, j]))
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
    
    for j = 1:Ny
        for i = 1:Nx
            Z[i, j] = sqrt(1 - X[i, j]^2 - Y[i, j]^2)
        end
    end
    
    time_end = now()
    time_elapsed = time_end - time_start
    
    if converged
        @printf("The SOR has converged after %d iterations.\n", iteration_final)
        @printf("The number of iterations required for convergence is %d.\n", iteration_final)
        @printf("Elapsed CPU time is %s.\n", time_elapsed)
    else
        @printf("The SOR has not converged after %d iterations.\n", iteration_final)
    end
    
    @printf("Average time per iteration is %.3g.\n", Dates.value(time_elapsed)/(1000 * iteration_final))
    
    if return_computational_coordinates
        return x, y, X, Y, Z
    else
        return X, Y, Z
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

    colors = [:blue, :red, :orange, :purple, :green, :darkcyan]
    alphas = [1, 1, 0.1125, 0.1125, 1, 0.1125]

    fig = Figure(resolution = (1500, 750))

    ax2D = Axis(fig[1, 1]; aspect = 1, title = "Cubed Sphere", axis_kwargs_2D...)
    ax3D = Axis3(fig[1, 2]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = "Cubed Sphere", 
                 axis_kwargs_3D...)
    
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
    colgap!(fig.layout, 40)
    
    save(filename, fig)
    
end


Nx, Ny = 32, 32

X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
filename = "conformal_cubed_sphere_2D_3D.png"
visualize_conformal_cubed_sphere_2D_3D(X, Y, Z, filename)

X, Y, Z = elliptic_quasi_conformal_cubed_sphere_coordinates(Nx, Ny)
filename = "elliptic_quasi_conformal_cubed_sphere_2D_3D.png"
visualize_conformal_cubed_sphere_2D_3D(X, Y, Z, filename)