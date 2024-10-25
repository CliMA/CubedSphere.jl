using Rotations
using CubedSphere
using DelimitedFiles
using CairoMakie


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


function visualize_conformal_cubed_sphere_panel_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, color)
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
    
    fig = Figure(resolution = (750, 750))
    
    ax2D = Axis(fig[1, 1]; aspect = 1, title = "Conformal Cubed Sphere Panel", axis_kwargs_2D...)
    hide_decorations && hidedecorations!(ax2D)
    wireframe!(ax2D, X, Y, Z, color = color)
    
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 40)

    return fig
end


function visualize_conformal_cubed_sphere_panel_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, color)
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
    
    fig = Figure(resolution = (750, 750))
    
    ax3D = Axis3(fig[1, 1]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)),
                 title = "Conformal Cubed Sphere Panel", axis_kwargs_3D...)
    hide_decorations && hidedecorations!(ax3D)
    wireframe!(ax3D, X, Y, Z, color = color)
    
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 40)

    return fig
end


function visualize_conformal_cubed_sphere_panel_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, color)
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
    
    fig = Figure(resolution = (1500, 750))
    
    ax2D = Axis(fig[1, 1]; aspect = 1, title = "Conformal Cubed Sphere Panel", axis_kwargs_2D...)
    ax3D = Axis3(fig[1, 2]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)),
                 title = "Conformal Cubed Sphere Panel", axis_kwargs_3D...)
    
    for ax in [ax2D, ax3D]
        hide_decorations && hidedecorations!(ax)
        wireframe!(ax, X, Y, Z, color = color)
    end
    
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 40)

    return fig
end


function visualize_conformal_cubed_sphere_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, colors;
                                             title = "Conformal Cubed Sphere: 2D Projection")
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

    fig = Figure(resolution = (750, 750))

    ax2D = Axis(fig[1, 1]; aspect = 1, title = title, axis_kwargs_2D...)
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
    
    return fig
end


function visualize_conformal_cubed_sphere_2D(X, Y, Z, filename; title = "Conformal Cubed Sphere: 2D Projection")
    axis_kwargs_2D = (xlabelsize = 22.5, ylabelsize = 22.5, xticklabelsize = 17.5, yticklabelsize = 17.5,
                      xticklabelpad = 10, yticklabelpad = 10, titlesize = 27.5, titlegap = 15, titlefont = :bold,
                      xlabel = "x", ylabel = "y")
    hide_decorations = false

    colors = [:orange, :red, :deepskyblue, :purple, :green, :blue]

    fig = Figure(resolution = (750, 750))

    ax2D = Axis(fig[1, 1]; aspect = 1, title = title, axis_kwargs_2D...)
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


function visualize_conformal_cubed_sphere_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, colors, alphas;
                                             title = "Conformal Cubed Sphere: 3D View")
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

    fig = Figure(resolution = (750, 750))

    ax3D = Axis3(fig[1, 1]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = title, axis_kwargs_3D...)
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
    
    return fig
end


function visualize_conformal_cubed_sphere_3D(X, Y, Z, filename; title = "Conformal Cubed Sphere: 3D View")
    axis_kwargs_3D = (xlabelsize = 22.5, ylabelsize = 22.5, zlabelsize = 22.5, xticklabelsize = 17.5,
                      yticklabelsize = 17.5, zticklabelsize = 17.5, xticklabelpad = 10, yticklabelpad = 10,
                      zticklabelpad = 10, titlesize = 27.5, titlegap = 15, titlefont = :bold, xlabel = "x",
                      ylabel = "y", zlabel = "z")
    hide_decorations = false

    colors = [:orange, :red, :deepskyblue, :purple, :green, :blue]
    alphas = [1, 1, 0.1125, 0.1125, 1, 0.1125]

    fig = Figure(resolution = (750, 750))

    ax3D = Axis3(fig[1, 1]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = title, axis_kwargs_3D...)
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


function visualize_conformal_cubed_sphere_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, colors,
                                                alphas; title_2D = "Conformal Cubed Sphere: 2D Projection",
                                                title_3D = "Conformal Cubed Sphere: 3D View")
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

    fig = Figure(resolution = (1500, 750))

    ax2D = Axis(fig[1, 1]; aspect = 1, title = title_2D, axis_kwargs_2D...)
    ax3D = Axis3(fig[1, 2]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = title_3D,
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
    
    return fig
end


function visualize_conformal_cubed_sphere_2D_3D(X, Y, Z, filename; title_2D = "Conformal Cubed Sphere: 2D Projection",
                                                title_3D = "Conformal Cubed Sphere: 3D View")
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

    fig = Figure(resolution = (1500, 750))

    ax2D = Axis(fig[1, 1]; aspect = 1, title = title_2D, axis_kwargs_2D...)
    ax3D = Axis3(fig[1, 2]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = title_3D,
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
    colgap!(fig.layout, 200)
    
    save(filename, fig)
end
