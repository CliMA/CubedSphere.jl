using CairoMakie
using CubedSphere
using Dates
using DelimitedFiles
using FilePaths
using Printf
using LinearAlgebra
using Random
using Rotations
using Statistics: mean, norm


function write_output_to_file_1D(output_directory, x, y, filename)
    cwd = pwd()
    path = joinpath(cwd, output_directory)
    if !isdir(path)
        mkdir(path)
    end
    cd(path)
    
    filename *= ".curve"
    outputfile = open(filename, "w")
    
    write(outputfile, "#phi\n")
    for i in eachindex(x)
        write(outputfile, string(x[i], " ", y[i], "\n"))
        # The above line is equivalent to println(outputfile, "$(x[i]) $(y[i])")
    end
    
    close(outputfile)
    cd(cwd)
end


function read_output_from_file_1D(output_directory, filename)
    cwd = pwd()
    path = joinpath(cwd, output_directory)
    if !isdir(path)
        mkdir(path)
    end
    cd(path)
    
    data = []
    count = 1
    open(filename, "r") do infile
        for line in eachline(infile)
            if count != 1
                push!(data, line)
            end
            count += 1
        end
    end
    data = readdlm(IOBuffer(join(data, "\n")))
    
    N = size(data, 1)
    x = zeros(N)
    y = zeros(N)
    for i in 1:N
        x[i] = data[i,1]
        y[i] = data[i,2]
    end
    
    cd(cwd)
    
    return (x, y)
end


function create_single_line_or_scatter_plot(resolution, plot_type, x, y, axis_kwargs, title, plot_kwargs, file_name;
                                            specify_x_limits = false, x_limits = [0, 0], specify_y_limits = false,
                                            y_limits = [0, 0], tight_x_axis = false, tight_y_axis = false,
                                            format = ".png")
    fig = Figure(resolution = resolution)
    ax = Axis(fig[1,1]; axis_kwargs...)

    if plot_type == "line_plot"
        lines!(ax, x, y, linewidth = plot_kwargs.linewidth, color = plot_kwargs.linecolor)
    elseif plot_type == "scatter_plot"
        scatter!(ax, x, y, marker = plot_kwargs.marker, markersize = plot_kwargs.markersize,
                 color = plot_kwargs.linecolor)
    elseif plot_type == "scatter_line_plot"
        scatterlines!(ax, x, y, linewidth = plot_kwargs.linewidth, marker = plot_kwargs.marker,
                      markersize = plot_kwargs.markersize, color = plot_kwargs.linecolor)
    end
    ax.title = title

    if specify_x_limits
        xlims!(ax, x_limits...)
    elseif tight_x_axis
        xlims!(ax, extrema(x)...)
    end

    if specify_y_limits
        ylims!(ax, y_limits...)
    elseif tight_y_axis
        ylims!(ax, extrema(y)...)
    end

    save(file_name * format, fig)
end


function test_create_single_line_or_scatter_plot()
    x = range(0, 2π, length = 100)
    y = sin.(x)

    resolution = (850, 750)

    axis_kwargs = (xlabel = "x", ylabel = "sin(x)", xlabelsize = 22.5, ylabelsize = 22.5, xticklabelsize = 17.5,
                   yticklabelsize = 17.5, xlabelpadding = 10, ylabelpadding = 10, aspect = 1.0, titlesize = 27.5,
                   titlegap = 15, titlefont = :bold)
    title = "sin(x) vs x"
    plot_kwargs = (linewidth = 2, linecolor = :black, marker = :rect, markersize = 10)

    plot_types = ["line_plot", "scatter_plot", "scatter_line_plot"]
    file_names = ["LinePlotExample", "ScatterPlotExample", "ScatterLinePlotExample"]

    for i in 1:3
        plot_type = plot_types[i]
        file_name = file_names[i]
        create_single_line_or_scatter_plot(resolution, plot_type, x, y, axis_kwargs, title, plot_kwargs, file_name)
    end
end


test_create_single_line_or_scatter_plot()


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
    deviation_from_isotropy = zeros(Nx-1, Ny-1)

    for i in 1:Nx-1
        for j in 1:Ny-1
            deviation_from_isotropy[i, j] = compute_deviation_from_isotropy_at_point(X, Y, i, j)
        end
    end

    return norm(deviation_from_isotropy)
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
    angle1 = abs(acos(clamp(dot(normalize(v1), normalize(-v4)), -1, 1)))
    angle3 = abs(acos(clamp(dot(normalize(v3), normalize(-v2)), -1, 1)))

    # Compute area using triangles formed by the cell.
    s1 = 0.5 * norm(v1) * norm(v4) * sin(angle1)
    s3 = 0.5 * norm(v3) * norm(v2) * sin(angle3)

    cell_area = s1 + s3

    return cell_area
end


function compute_cell_areas(X, Y)
    Nx, Ny = size(X)
    cell_areas = zeros(Nx-1, Ny-1)

    for i in 1:Nx-1
        for j in 1:Ny-1
            cell_areas[i, j] = compute_cell_area(X, Y, i, j)
        end
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


function visualize_cubed_sphere_2D(X, Y, Z, filename)
    axis_kwargs_2D = (xlabelsize = 22.5, ylabelsize = 22.5, xticklabelsize = 17.5, yticklabelsize = 17.5, 
                      xticklabelpad = 10, yticklabelpad = 10, titlesize = 27.5, titlegap = 15, titlefont = :bold, 
                      xlabel = "x", ylabel = "y")
    hide_decorations = false
    
    colors = [:orange, :red, :deepskyblue, :purple, :green, :blue]

    fig = Figure(resolution = (750, 750))

    # ax2D = Axis(fig[1, 1]; aspect = 1, title = "Cubed Sphere", axis_kwargs_2D...)
    ax2D = Axis(fig[1, 1]; aspect = 1, title = "", axis_kwargs_2D...)
    hide_decorations && hidedecorations!(ax2D)
    wireframe!(ax2D, X, Y, Z, color = colors[1])

    rotations = (RotX(π/2), RotX(-π/2), RotY(π/2), RotY(-π/2), RotX(π))
    
    for (i, R) in enumerate(rotations)
        X′ = similar(X)
        Y′ = similar(Y)
        Z′ = similar(Z)
    
        for I in CartesianIndices(X)
            X′[I], Y′[I], Z′[I] = R * [X[I], Y[I], Z[I]]
        end
    
        wireframe!(ax2D, X′, Y′, Z′, color = colors[i + 1])
    end
    
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 40)
    
    save(filename, fig)
end


function visualize_cubed_sphere_3D(X, Y, Z, filename)
    axis_kwargs_3D = (xlabelsize = 22.5, ylabelsize = 22.5, zlabelsize = 22.5, xticklabelsize = 17.5, 
                      yticklabelsize = 17.5, zticklabelsize = 17.5, xticklabelpad = 10, yticklabelpad = 10, 
                      zticklabelpad = 10, titlesize = 27.5, titlegap = 15, titlefont = :bold, xlabel = "x", 
                      ylabel = "y", zlabel = "z")
    hide_decorations = false

    colors = [:orange, :red, :deepskyblue, :purple, :green, :blue]
    alphas = [1, 1, 0.1125, 0.1125, 1, 0.1125]
    
    fig = Figure(resolution = (800, 750))

    ax3D = Axis3(fig[1, 1]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = "Cubed Sphere", 
                 axis_kwargs_3D...)
    hide_decorations && hidedecorations!(ax3D)
    
    wireframe!(ax3D, X, Y, Z, color = colors[1], alpha = alphas[1])
    
    rotations = (RotX(π/2), RotX(-π/2), RotY(π/2), RotY(-π/2), RotX(π))
    
    for (i, R) in enumerate(rotations)
        X′ = similar(X)
        Y′ = similar(Y)
        Z′ = similar(Z)
    
        for I in CartesianIndices(X)
            X′[I], Y′[I], Z′[I] = R * [X[I], Y[I], Z[I]]
        end
    
        wireframe!(ax3D, X′, Y′, Z′, color = colors[i + 1], alpha = alphas[i + 1])
    end
    
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 40)
    
    save(filename, fig)
end


function visualize_cubed_sphere_2D_3D(X, Y, Z, filename)
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


function compute_model_diagnostics(X, Y, minimum_reference_cell_area)
    cell_areas = compute_cell_areas(X, Y)
    normalized_minimum_cell_width = sqrt(minimum(cell_areas)/minimum_reference_cell_area)
    
    deviation_from_isotropy = compute_deviation_from_isotropy(X, Y)

    model_diagnostics = vcat(normalized_minimum_cell_width, deviation_from_isotropy)
    
    return model_diagnostics
end


function compute_weighted_model_diagnostics(Nx, Ny, model_diagnostics)
    normalized_minimum_cell_width = model_diagnostics[1]
    deviation_from_isotropy = model_diagnostics[2]

    weights = specify_weights_for_model_diagnostics()

    weighted_model_diagnostics = vcat(weights[1] * normalized_minimum_cell_width,
                                      weights[2] * deviation_from_isotropy)

    return weighted_model_diagnostics
end


function forward_map(Nx, Ny, spacing_type, θ)
    θ_limits = specify_parameter_limits(spacing_type)

    for i in 1:lastindex(θ)
        θ[i] = clamp(θ[i], θ_limits[i][1], θ_limits[i][2])
    end
    
    x_reference, y_reference, X_reference, Y_reference, Z_reference = conformal_cubed_sphere_coordinates(Nx, Ny)
    cell_areas = compute_cell_areas(X_reference, Y_reference)
    minimum_reference_cell_area = minimum(cell_areas)
    
    x, y, X, Y, Z = (
    conformal_cubed_sphere_coordinates(Nx, Ny; non_uniform_spacing = true, spacing_type = spacing_type,
                                       ratio_raised_to_Nx_minus_one = θ[1], k₀ByNx = θ[1]))
    
    model_diagnostics = compute_model_diagnostics(X, Y, minimum_reference_cell_area)

    weighted_model_diagnostics = compute_weighted_model_diagnostics(Nx, Ny, model_diagnostics)

	return weighted_model_diagnostics
end


function specify_ideal_weighted_model_diagnostics(Nx, Ny)
    normalized_minimum_cell_width = 4

    deviation_from_isotropy = 0.0

    weights = specify_weights_for_model_diagnostics()

    ideal_weighted_model_diagnostics = vcat(weights[1] * normalized_minimum_cell_width,
                                            weights[2] * deviation_from_isotropy)

    return ideal_weighted_model_diagnostics
end


function optimize!(Nx, Ny, spacing_type, θ; nIterations = 10, Δt = 1)
    ideal_data = specify_ideal_weighted_model_diagnostics(Nx, Ny)
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

		# Evaluating the forward map for all ensemble members. This is the most expensive step because it need to  run
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


function optimized_non_uniform_conformal_cubed_sphere(Nx, Ny, spacing_type)
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


function visualize_optimized_non_uniform_conformal_cubed_sphere(Nx, Ny, spacing_type, optimized)
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
    visualize_cubed_sphere_2D(X, Y, Z, "conformal_cubed_sphere_2D.png")
    visualize_cubed_sphere_3D(X, Y, Z, "conformal_cubed_sphere_3D.png")
    visualize_cubed_sphere_2D_3D(X, Y, Z, "conformal_cubed_sphere_2D_3D.png")
    cell_areas = compute_cell_areas(X, Y)
    reference_minimum_cell_area = minimum(cell_areas)
    
    if optimized
        x, y, X, Y, Z = optimized_non_uniform_conformal_cubed_sphere(Nx, Ny, spacing_type)
    else
        x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny; non_uniform_spacing = true,
                                                           spacing_type = spacing_type)
    end
    
    if optimized
        file_name_suffix = "_optimized"
    else
        file_name_suffix = "_unoptimized"
    end
    visualize_cubed_sphere_2D(X, Y, Z,
                              "non_uniform_conformal_cubed_sphere_2D_" * spacing_type * file_name_suffix * ".png")
    visualize_cubed_sphere_3D(X, Y, Z,
                              "non_uniform_conformal_cubed_sphere_3D_" * spacing_type * file_name_suffix * ".png")
    visualize_cubed_sphere_2D_3D(X, Y, Z,
                                 "non_uniform_conformal_cubed_sphere_2D_3D_" * spacing_type * file_name_suffix * ".png")
    
    cell_areas = compute_cell_areas(X, Y)
    minimum_cell_area = minimum(cell_areas)
    print("The normalized minimum cell width of the non-uniform conformal cubed sphere for spacing = $spacing_type " 
          * "and optimized = $optimized is $(sqrt(minimum_cell_area/reference_minimum_cell_area))\n")
end 


for optimized in [false, true]
    for spacing_type in ["geometric", "exponential"]
        @info "Visualizing non-uniform conformal cubed sphere for spacing = $spacing_type and optimized = $optimized"
        N = 16
        Nx, Ny = N + 1, N + 1
        visualize_optimized_non_uniform_conformal_cubed_sphere(Nx, Ny, spacing_type, optimized)
    end
end


function minimum_cell_width_variation_with_resolution(spacing_type, optimized)
    resolutions = 100:50:1000
    normalized_minimum_cell_widths = zeros(length(resolutions))
    
    minimum_reference_cell_area = 0

    for (i, resolution) in enumerate(resolutions)
        Nx, Ny = resolution + 1, resolution + 1
        x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
        cell_areas = compute_cell_areas(X, Y)
        minimum_reference_cell_area = minimum(cell_areas)
        if optimized
            x, y, X, Y, Z = optimized_non_uniform_conformal_cubed_sphere(Nx, Ny, spacing_type)
        else
            x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny; non_uniform_spacing = true,
                                                               spacing_type = spacing_type)
        end
        cell_areas = compute_cell_areas(X, Y)
        minimum_cell_area = minimum(cell_areas)
        normalized_minimum_cell_widths[i] = sqrt(minimum_cell_area/minimum_reference_cell_area)
    end

    plot_resolution = (800, 750)
    plot_type = "scatter_line_plot"
    xlabel = "Number of cells in each direction of a cubed sphere panel"
    ylabel = "Normalized minimum cell width"
    axis_kwargs = (xlabel = "Number of cells in each direction of a cubed sphere panel",
                   ylabel = "Normalized minimum cell width", xlabelsize = 22.5, ylabelsize = 22.5,
                   xticklabelsize = 17.5, yticklabelsize = 17.5, xlabelpadding = 10, ylabelpadding = 10, aspect = 1.0,
                   titlesize = 27.5, titlegap = 15, titlefont = :bold)
    title = "Normalized Minimum Cell Width versus Resolution"
    plot_kwargs = (linewidth = 2, linecolor = :black, marker = :rect, markersize = 15)
    if spacing_type == "geometric"
        file_name_suffix_1 = "_Geometric"
    elseif spacing_type == "exponential"
        file_name_suffix_1 = "_Exponential"
    end
    if optimized
        file_name_suffix_2 = "_Optimized"
    else
        file_name_suffix_2 = "_Unoptimized"
    end
    file_name = "MinimumCellWidthVersusResolution" * file_name_suffix_1 * file_name_suffix_2
    
    write_output_to_file_1D(".", resolutions, normalized_minimum_cell_widths, file_name)
    resolutions, normalized_minimum_cell_widths = read_output_from_file_1D(".", file_name * ".curve")
    create_single_line_or_scatter_plot(plot_resolution, plot_type, resolutions[2:end],
                                       normalized_minimum_cell_widths[2:end], axis_kwargs, title, plot_kwargs,
                                       file_name)
end


for optimized in [false, true]
    for spacing_type in ["geometric", "exponential"]
        @info "Computing minimum cell width variation with resolution for spacing = $spacing_type " *
        "and optimized = $optimized"
        minimum_cell_width_variation_with_resolution(spacing_type, optimized)
    end
end
