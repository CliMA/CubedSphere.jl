# CubedSphere.jl Documentation

## Conformal cubed sphere mapping

Imagine a cube inscribed into a sphere. Using [`conformal_cubed_sphere_mapping`](@ref) we can map the face of a cube onto the sphere. [`conformal_cubed_sphere_mapping`](@ref) maps the face that corresponds to the sphere
sector that includes the North Pole, that is, the face of the cube is oriented normal
to the ``z`` axis and it is parametrized in the range ``(x, y) \in [-1, 1] \times [-1, 1]``.

We can visualize how this mapping looks like.

```@setup 1
using Rotations
using CairoMakie
CairoMakie.activate!(type = "svg")
using CubedSphere
```

```@example 1
using CairoMakie
using CubedSphere

N = 16

x = range(-1, 1, length=N)
y = range(-1, 1, length=N)

X = zeros(length(x), length(y))
Y = zeros(length(x), length(y))
Z = zeros(length(x), length(y))

for (j, y′) in enumerate(y), (i, x′) in enumerate(x)
    X[i, j], Y[i, j], Z[i, j] = conformal_cubed_sphere_mapping(x′, y′)
end

fig = Figure()

ax1 = Axis(fig[1, 1], aspect = 1)
ax2 = Axis3(fig[1, 2], aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)))

for ax in [ax1, ax2]
    hidedecorations!(ax)
    wireframe!(ax, X, Y, Z)
end

colsize!(fig.layout, 1, Auto(0.8))
colgap!(fig.layout, 40)

current_figure()
```

Above, we plotted both in a 2D projection and in 3D space.

We can then use [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl) to rotate the face of the sphere
that includes the North Pole and this way obtain all six faces of the sphere.

```@example 1
using Rotations

fig = Figure()

ax1 = Axis(fig[1, 1], aspect=1)
ax2 = Axis3(fig[1, 2], aspect=(1, 1, 1), limits=((-1, 1), (-1, 1), (-1, 1)))

for ax in [ax1, ax2]
    wireframe!(ax, X, Y, Z)
end

rotations = (RotX(π/2), RotX(-π/2), RotY(π/2), RotY(-π/2), RotX(π))

for R in rotations
    X′ = similar(X)
    Y′ = similar(Y)
    Z′ = similar(Z)

    for I in CartesianIndices(X)
        X′[I], Y′[I], Z′[I] = R * [X[I], Y[I], Z[I]]
    end

    for ax in [ax1, ax2]
        wireframe!(ax, X′, Y′, Z′)
    end
end

colsize!(fig.layout, 1, Auto(0.8))
colgap!(fig.layout, 40)

current_figure()
```
