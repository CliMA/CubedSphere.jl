# CubedSphere.jl Documentation

## Conformal cubed sphere mapping

The conformal method for projecting the cube on the sphere was first described the paper by [Rancic-etal-1996](@citet).

> Rančić et al., (1996). A global shallow-water model using an expanded spherical cube - Gnomonic versus conformal coordinates, _Quarterly Journal of the Royal Meteorological Society_, **122 (532)**, 959-982. doi:[10.1002/qj.49712253209](https://doi.org/10.1002/qj.49712253209)

Imagine a cube inscribed into a sphere. Using [`conformal_cubed_sphere_mapping`](@ref) we can map the face of the
cube onto the sphere. [`conformal_cubed_sphere_mapping`](@ref) maps the face that corresponds to the sphere's
sector that includes the North Pole, that is, the face of the cube is oriented normal to the ``z`` axis. This cube's
face is parametrized with orthogonal coordinates ``(x, y) \in [-1, 1] \times [-1, 1]`` with ``(x, y) = (0, 0)`` being
in the center of the cube's face, that is on the ``z`` axis.

We can visualize the mapping.

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

axis_kwargs_2D = (xlabelsize = 22.5, ylabelsize = 22.5, xticklabelsize = 17.5, yticklabelsize = 17.5, 
                  xticklabelpad = 10, yticklabelpad = 10, titlesize = 27.5, titlegap = 15, titlefont = :bold, 
                  xlabel = "x", ylabel = "y")

axis_kwargs_3D = (xlabelsize = 22.5, ylabelsize = 22.5, zlabelsize = 22.5, xticklabelsize = 17.5, yticklabelsize = 17.5,
                  zticklabelsize = 17.5, xticklabelpad = 10, yticklabelpad = 10, zticklabelpad = 10, titlesize = 27.5,
                  titlegap = 15, titlefont = :bold, xlabel = "x", ylabel = "y", zlabel = "z")

color = :blue
    
fig = Figure(resolution = (1500, 750))

ax2D = Axis(fig[1, 1]; aspect = 1, title = "Cubed Sphere Panel", axis_kwargs_2D...)
ax3D = Axis3(fig[1, 2]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = "Cubed Sphere Panel", 
             axis_kwargs_3D...)

for ax in [ax2D, ax3D]
    wireframe!(ax, X, Y, Z, color = color)
end

colsize!(fig.layout, 1, Auto(0.8))
colgap!(fig.layout, 40)

current_figure()
```

Above, we plotted the mapping from the cube's face onto the sphere both in a 2D projection (e.g., overlooking the sphere down to its North Pole) and in 3D space.

We can then use [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl) to rotate the face of the sphere that includes the North Pole and this way obtain all six faces of the sphere.

```@example 1
using Rotations

colors = [:purple, :red, :orange, :cyan, :green, :blue]
alphas = [1, 1, 0.1125, 0.1125, 1, 0.1125]

fig = Figure(resolution = (1500, 750))

ax2D = Axis(fig[1, 1]; aspect = 1, title = "Cubed Sphere", axis_kwargs_2D...)
ax3D = Axis3(fig[1, 2]; aspect = (1, 1, 1), limits = ((-1, 1), (-1, 1), (-1, 1)), title = "Cubed Sphere", 
             axis_kwargs_3D...)

for ax in [ax2D, ax3D]
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

current_figure()
```
