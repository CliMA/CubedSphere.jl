using Rotations
using GLMakie
using CubedSphere

x = range(-1, 1, length=20)
y = range(-1, 1, length=20)

X = zeros(length(x), length(y))
Y = zeros(length(x), length(y))
Z = zeros(length(x), length(y))

for (i, x′) in enumerate(x), (j, y′) in enumerate(y)
    X[i, j], Y[i, j], Z[i, j] = CubedSphere.conformal_cubed_sphere_mapping(x′, y′)
end

fig = wireframe(X, Y, Z, show_axis=false)

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
