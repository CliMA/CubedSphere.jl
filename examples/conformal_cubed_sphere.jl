include("file_io_and_plotting.jl")


function single_line_or_scatter_plot_example()
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


single_line_or_scatter_plot_example()


function conformal_cubed_sphere_2D_3D_visualization_example()
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
    colors = [:orange, :red, :deepskyblue, :purple, :green, :blue]
    alphas = [1, 1, 0.1125, 0.1125, 1, 0.1125]

    fig = visualize_conformal_cubed_sphere_panel_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, color)
    save("conformal_cubed_sphere_panel_2D.png", fig)

    fig = visualize_conformal_cubed_sphere_panel_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, color)
    save("conformal_cubed_sphere_panel_3D.png", fig)

    fig = visualize_conformal_cubed_sphere_panel_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, color)
    save("conformal_cubed_sphere_panel_2D_3D.png", fig)

    fig = visualize_conformal_cubed_sphere_2D(Nx, Ny, axis_kwargs_2D, hide_decorations, colors)
    save("conformal_cubed_sphere_2D.png", fig)

    fig = visualize_conformal_cubed_sphere_3D(Nx, Ny, axis_kwargs_3D, hide_decorations, colors, alphas)
    save("conformal_cubed_sphere_3D.png", fig)

    fig = visualize_conformal_cubed_sphere_2D_3D(Nx, Ny, axis_kwargs_2D, axis_kwargs_3D, hide_decorations, colors,
                                                 alphas)
    save("conformal_cubed_sphere_2D_3D.png", fig)
end


conformal_cubed_sphere_2D_3D_visualization_example()


function visualize_optimized_non_uniform_conformal_cubed_sphere(Nx, Ny; spacing_type = "geometric", optimized = false)
    x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
    visualize_conformal_cubed_sphere_2D(X, Y, Z, "conformal_cubed_sphere_2D.png")
    visualize_conformal_cubed_sphere_3D(X, Y, Z, "conformal_cubed_sphere_3D.png")
    visualize_conformal_cubed_sphere_2D_3D(X, Y, Z, "conformal_cubed_sphere_2D_3D.png")
    cell_areas = compute_cell_areas(X, Y, Z)
    reference_minimum_cell_area = minimum(cell_areas)
    
    if optimized
        x, y, X, Y, Z = optimized_non_uniform_conformal_cubed_sphere_coordinates(Nx, Ny, spacing_type)
    else
        x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny; non_uniform_spacing = true,
                                                           spacing_type = spacing_type)
    end
    
    if spacing_type == "geometric"
        spacing_type_title = "Geometric"
    elseif spacing_type == "exponential"
        spacing_type_title = "Exponential"
    end

    if optimized
        file_name_suffix = "_optimized"
        title_2D = "CCS with Adaptive $spacing_type_title Spacing and\nOptimized with EKI: 2D Projection"
        title_3D = "CCS with Adaptive $spacing_type_title Spacing\nand Optimized with EKI: 3D View"
    else
        file_name_suffix = "_unoptimized"
        title_2D = "CCS with Adaptive $spacing_type_title Spacing:\n2D Projection"
        title_3D = "CCS with Adaptive $spacing_type_title Spacing:\n3D View"
    end
    visualize_conformal_cubed_sphere_2D(
    X, Y, Z, "non_uniform_conformal_cubed_sphere_2D_" * spacing_type * file_name_suffix * ".png"; title = title_2D)
    visualize_conformal_cubed_sphere_3D(
    X, Y, Z, "non_uniform_conformal_cubed_sphere_3D_" * spacing_type * file_name_suffix * ".png"; title = title_3D)
    visualize_conformal_cubed_sphere_2D_3D(
    X, Y, Z, "non_uniform_conformal_cubed_sphere_2D_3D_" * spacing_type * file_name_suffix * ".png";
    title_2D = title_2D, title_3D = title_3D)
    
    cell_areas = compute_cell_areas(X, Y, Z)
    minimum_cell_area = minimum(cell_areas)
    print("The normalized minimum cell width of the non-uniform conformal cubed sphere for spacing = $spacing_type " 
          * "and optimized = $optimized is $(sqrt(minimum_cell_area/reference_minimum_cell_area))\n")
end 


for optimized in [false, true]
    for spacing_type in ["geometric", "exponential"]
        @info "Visualizing non-uniform conformal cubed sphere for spacing = $spacing_type and optimized = $optimized"
        N = 16
        Nx, Ny = N + 1, N + 1
        visualize_optimized_non_uniform_conformal_cubed_sphere(Nx, Ny; spacing_type = spacing_type, optimized=optimized)
    end
end


function minimum_cell_width_variation_with_resolution(spacing_type, optimized)
    resolutions = 100:50:1000
    normalized_minimum_cell_widths = zeros(length(resolutions))
    
    minimum_reference_cell_area = 0

    for (i, resolution) in enumerate(resolutions)
        Nx, Ny = resolution + 1, resolution + 1
        x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny)
        cell_areas = compute_cell_areas(X, Y, Z)
        minimum_reference_cell_area = minimum(cell_areas)
        if optimized
            x, y, X, Y, Z = optimized_non_uniform_conformal_cubed_sphere_coordinates(Nx, Ny, spacing_type)
        else
            x, y, X, Y, Z = conformal_cubed_sphere_coordinates(Nx, Ny; non_uniform_spacing = true,
                                                               spacing_type = spacing_type)
        end
        cell_areas = compute_cell_areas(X, Y, Z)
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


compute_and_plot_minimum_cell_width_variation_with_resolution = false

if compute_and_plot_minimum_cell_width_variation_with_resolution
    for optimized in [false, true]
        for spacing_type in ["geometric", "exponential"]
            @info "Computing minimum cell width variation with resolution for spacing = $spacing_type " *
            "and optimized = $optimized"
            minimum_cell_width_variation_with_resolution(spacing_type, optimized)
        end
    end
end
