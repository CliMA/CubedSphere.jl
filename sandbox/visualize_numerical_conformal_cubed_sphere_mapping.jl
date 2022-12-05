using Rotations
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

N = 16

x = range(-1, 1, length=N)
y = range(-1, 1, length=N)

X = zeros(length(x), length(y))
Y = zeros(length(x), length(y))
Z = zeros(length(x), length(y))

for (i, x′) in enumerate(x), (j, y′) in enumerate(y)
    X[i, j], Y[i, j], Z[i, j] = conformal_cubed_sphere_mapping(x′, y′)
end

fig = Figure()
ax = Axis(fig[1, 1], aspect=1)

wireframe!(ax, X, Y, Z, show_axis=false)

rotations = (RotX(π/2), RotX(-π/2), RotY(π/2), RotY(-π/2), RotX(π))

for R in rotations
    X′ = similar(X)
    Y′ = similar(Y)
    Z′ = similar(Z)

    for I in CartesianIndices(X)
        X′[I], Y′[I], Z′[I] = R * [X[I], Y[I], Z[I]]
    end

    wireframe!(X′, Y′, Z′)
end

display(fig)
