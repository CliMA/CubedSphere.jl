using CubedSphere
using GLMakie

x = range(-1, 1, length=51)
y = range(-1, 1, length=51)
X = zeros(length(x), length(y))
Y = zeros(length(x), length(y))
Z = zeros(length(x), length(y))

for (i, x′) in enumerate(x), (j, y′) in enumerate(y)
    X[i, j], Y[i, j], Z[i, j] = CubedSphere.conformal_cubed_sphere_mapping(x′, y′)
end

wireframe(X, Y, Z)

