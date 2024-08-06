using Rotations
using CairoMakie
using CubedSphere

function conformal_cubed_sphere_coordinates(Nx, Ny)

    x = range(-1, 1, length = Nx)
    y = range(-1, 1, length = Ny)
    
    X = zeros(length(x), length(y))
    Y = zeros(length(x), length(y))
    Z = zeros(length(x), length(y))
    
    for (j, y′) in enumerate(y), (i, x′) in enumerate(x)
        X[i, j], Y[i, j], Z[i, j] = conformal_cubed_sphere_mapping(x′, y′)
    end
    
    return X, Y, Z
    
end

function visualize_cubed_sphere_panel_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, color)

    X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
    
    fig = Figure(resolution = (750, 750))
    
    ax2D = Axis(fig[1, 1]; aspect = 1, title = "Cubed Sphere Panel", axis_kwargs_2D...)
    hide_decorations && hidedecorations!(ax2D)
    wireframe!(ax2D, X, Y, Z, color = color)
    
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 40)

    return fig

end

function visualize_cubed_sphere_panel_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, color)

    X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
    
    fig = Figure(resolution = (750, 750))
    
    ax3D = Axis3(fig[1, 1]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = "Cubed Sphere Panel", 
                 axis_kwargs_3D...)
    hide_decorations && hidedecorations!(ax3D)
    wireframe!(ax3D, X, Y, Z, color = color)
    
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 40)

    return fig

end

function visualize_cubed_sphere_panel_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, color)

    X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
    
    fig = Figure(resolution = (1500, 750))
    
    ax2D = Axis(fig[1, 1]; aspect = 1, title = "Cubed Sphere Panel", axis_kwargs_2D...)
    ax3D = Axis3(fig[1, 2]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = "Cubed Sphere Panel", 
                 axis_kwargs_3D...)
    
    for ax in [ax2D, ax3D]
        hide_decorations && hidedecorations!(ax)
        wireframe!(ax, X, Y, Z, color = color)
    end
    
    colsize!(fig.layout, 1, Auto(0.8))
    colgap!(fig.layout, 40)

    return fig

end

function visualize_cubed_sphere_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, colors)
    
    X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

    fig = Figure(resolution = (750, 750))

    ax2D = Axis(fig[1, 1]; aspect = 1, title = "Cubed Sphere", axis_kwargs_2D...)
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

function visualize_cubed_sphere_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, colors, alphas)
    
    X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

    fig = Figure(resolution = (750, 750))

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
    
    return fig

end

function visualize_cubed_sphere_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, colors, alphas)
    
    X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)

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
    
    return fig

end

function test_visualize_cubed_sphere_2D_3D()

    Nx, Ny = 16, 16
    axis_kwargs_2D = (xlabelsize = 22.5, ylabelsize = 22.5, xticklabelsize = 17.5, yticklabelsize = 17.5,
                      xticklabelpad = 10, yticklabelpad = 10, titlesize = 27.5, titlegap = 15, titlefont = :bold,
                      xlabel = "x", ylabel = "y")
    axis_kwargs_3D = (xlabelsize = 22.5, ylabelsize = 22.5, zlabelsize = 22.5, xticklabelsize = 17.5,
                      yticklabelsize = 17.5, zticklabelsize = 17.5, xticklabelpad = 10, yticklabelpad = 10,
                      zticklabelpad = 10, titlesize = 27.5, titlegap = 15, titlefont = :bold, xlabel = "x",
                      ylabel = "y", zlabel = "z")
    hide_decorations = false
    color = :blue
    colors = [:purple, :red, :orange, :cyan, :green, :blue]
    alphas = [1, 1, 0.1125, 0.1125, 1, 0.1125]

    fig = visualize_cubed_sphere_panel_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, color)
    save("cubed_sphere_panel_2D.png", fig)

    fig = visualize_cubed_sphere_panel_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, color)
    save("cubed_sphere_panel_3D.png", fig)

    fig = visualize_cubed_sphere_panel_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, color)
    save("cubed_sphere_panel_2D_3D.png", fig)

    fig = visualize_cubed_sphere_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, colors)
    save("cubed_sphere_2D.png", fig)

    fig = visualize_cubed_sphere_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, colors, alphas)
    save("cubed_sphere_3D.png", fig)

    fig = visualize_cubed_sphere_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, colors, alphas)
    save("cubed_sphere_2D_3D.png", fig)

end

test_visualize_cubed_sphere_2D_3D()