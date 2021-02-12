using CairoMakie

"""
    visualize_conformal_mapping(w; x_min, y_min, x_max, y_max, n_lines, n_samples, filepath="mapping.png")

Visualize a complex mapping w(z) from the rectangle `x ∈ [x_min, x_max]`, `y ∈ [y_min, y_max]` using `n_lines`
equally spaced lines in x and y each evaluated at `n_samples` sample points.
"""
function visualize_conformal_mapping(w; x_min, y_min, x_max, y_max, n_lines, n_samples, filepath="mapping.png")
    z₁ = [[x + im * y for x in range(x_min, x_max, length=n_samples)] for y in range(y_min, y_max, length=n_lines)]
    z₂ = [[x + im * y for y in range(y_min, y_max, length=n_samples)] for x in range(x_min, x_max, length=n_lines)]
    zs = cat(z₁, z₂, dims=1)
    ws = [w.(z) for z in zs]

    fig = Figure(resolution=(1920, 1080))
    ax = fig[1, 1] = Axis(fig, xlabel="Re{w}", ylabel="Im{w}")
    [lines!(ax, real(z), imag(z), color=RGBAf0(0, 0, 0, 0.25)) for z in zs]
    [lines!(ax, real(w), imag(w), color="dodgerblue2") for w in ws]
    ax.aspect = AxisAspect(1)
    save(filepath, fig)
end
