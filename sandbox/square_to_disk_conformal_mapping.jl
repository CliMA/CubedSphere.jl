using GLMakie
using CubedSphere

"""
    visualize_conformal_mapping(w; x_min, y_min, x_max, y_max, n_lines, n_samples, filepath="mapping.png")

Visualize a complex mapping ``w(z)`` from the rectangle `x ∈ [x_min, x_max]`, `y ∈ [y_min, y_max]` using `n_lines`
equally spaced lines in ``x`` and ``y`` each evaluated at `n_samples` sample points.
"""
function visualize_conformal_mapping(w; x_min, y_min, x_max, y_max, n_lines, n_samples, filepath="mapping.png")

    z₁ = [[x + im * y for x in range(x_min, x_max, length=n_samples)] for y in range(y_min, y_max, length=n_lines)]
    z₂ = [[x + im * y for y in range(y_min, y_max, length=n_samples)] for x in range(x_min, x_max, length=n_lines)]

    zs = cat(z₁, z₂, dims=1)
    ws = [w.(z) for z in zs]

    fig = Figure(resolution=(1200, 1200))
    ax = Axis(fig[1, 1];
              xlabel = "Re{w}",
              ylabel = "Im{w}",
              aspect = 1)

    [lines!(ax, real(z), imag(z), color = :grey) for z in zs]
    [lines!(ax, real(w), imag(w), color = :blue) for w in ws]

    display(fig)

    save(filepath, fig)
end

w(z) = (1 - cn(z, 1/2)) / sn(z, 1/2)

visualize_conformal_mapping(w;
                            x_min = -1.5,
                            y_min = -1.5,
                            x_max =  1.5,
                            y_max =  1.5,
                            n_lines = 31,
                            n_samples = 101,
                            filepath = "square_to_disk_conformal_mapping.png")
