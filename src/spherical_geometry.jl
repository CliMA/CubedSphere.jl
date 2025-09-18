using Distances
using LinearAlgebra

"""
    cartesian_to_lat_lon(x, y, z)
    cartesian_to_lat_lon(X)

Convert 3D Cartesian coordinates `(x, y, z)` or a 3-element Cartesian vector `X = (x, y, z)` on the sphere to
latitude–longitude. Returns a tuple `(latitude, longitude)` in degrees.

- Latitude is the angle measured from the equatorial plane (`z = 0`).
- Longitude is measured anti-clockwise (eastward) from the `x`-axis (`y = 0`) about the `z`-axis.

# Arguments
- `x, y, z`: Cartesian coordinates (numbers), **or**
- `X`: 3-element Cartesian vector.

# Returns
- `(latitude, longitude)`: Latitude and longitude angles in degrees.

# Examples
Find latitude–longitude of the North Pole:

```jldoctest 1
julia> using CubedSphere

julia> x, y, z = (0, 0, 6.4e6); # Cartesian coordinates of the North Pole [in meters]

julia> cartesian_to_lat_lon(x, y, z)
(90.0, 0.0)
```
Let's confirm that for few points on the unit sphere we get the answers we expect.

```jldoctest 1
julia> cartesian_to_lat_lon(√2/4, -√2/4, √3/2)
(59.99999999999999, -45.0)

julia> cartesian_to_lat_lon(-√6/4, √2/4, -√2/2)
(-45.00000000000001, 150.0)
```
"""
cartesian_to_lat_lon(x, y, z) = cartesian_to_latitude(x, y, z), cartesian_to_longitude(x, y, z)

function cartesian_to_lat_lon(X)
    x, y, z = X
    return cartesian_to_lat_lon(x, y, z)
end

"""
    cartesian_to_latitude(x, y, z)

Convert Cartesian coordinates `(x, y, z)` to latitude (in degrees) on the sphere.
"""
cartesian_to_latitude(x, y, z) = atand(z, hypot(x, y))

"""
    cartesian_to_longitude(x, y, z)

Convert Cartesian coordinates `(x, y, z)` to longitude (in degrees) on the sphere.
"""
cartesian_to_longitude(x, y, z) = atand(y, x)

"""
    lat_lon_to_cartesian(φ, λ; radius = 1)

Convert `(latitude, longitude)` coordinates (in degrees) to Cartesian coordinates `(x, y, z)` on the sphere.

# Arguments
- `φ`: Latitude in degrees.
- `λ`: Longitude in degrees.
- `radius`: Sphere radius (optional). Default is `1`.

# Returns
- `(x, y, z)`: Cartesian coordinates on the sphere.

# Examples
Find the Cartesian coordinates of the North Pole on a unit sphere:

```jldoctest 1
julia> using CubedSphere

julia> lat_lon_to_cartesian(90, 0)
(0.0, 0.0, 1.0)
```
Find the Cartesian coordinates of a point on the equator with longitude 90°E:

```jldoctest 1
julia> lat_lon_to_cartesian(0, 90)
(0.0, 1.0, 0.0)
```
"""
function lat_lon_to_cartesian(φ, λ; radius = 1)
    abs(φ) > 90 && error("Latitude φ must be within -90 ≤ φ ≤ 90 degrees.")
    return (lat_lon_to_x(φ, λ; radius), lat_lon_to_y(φ, λ; radius), lat_lon_to_z(φ; radius))
end

"""
lat_lon_to_x(φ, λ; radius = 1)

Convert (latitude, longitude) coordinates (in degrees) to Cartesian coordinate x on the sphere.
"""
lat_lon_to_x(φ, λ; radius = 1) = radius * cosd(λ) * cosd(φ)

"""
lat_lon_to_y(φ, λ; radius = 1)

Convert (latitude, longitude) coordinates (in degrees) to Cartesian coordinate y on the sphere.
"""
lat_lon_to_y(φ, λ; radius = 1) = radius * sind(λ) * cosd(φ)

"""
lat_lon_to_z(φ; radius = 1)

Convert (latitude, longitude) coordinates (in degrees) to Cartesian coordinate z on the sphere.
"""
lat_lon_to_z(φ; radius = 1) = radius * sind(φ)

"""
    turning_angle(λ₁, φ₁, λ₂, φ₂)

Compute the **signed** turning angle (in degrees) between the unit tangent vectors at the endpoints of the great-circle
arc connecting two points `(λ₁, φ₁)` and `(λ₂, φ₂)` on the unit sphere.

# Arguments
- `λ₁, φ₁`: Longitude and latitude of the first point (in degrees).
- `λ₂, φ₂`: Longitude and latitude of the second point (in degrees).

# Returns
- Signed turning angle (in degrees) in `(-180, 180]`.

# Notes
- A positive angle corresponds to a counter-clockwise rotation from the tangent at `(λ₁, φ₁)` to the tangent at
  `(λ₂, φ₂)` about the great-circle normal.
- The result is undefined for coincident or antipodal points (an error is thrown).
"""

function turning_angle(φ₁, λ₁, φ₂, λ₂)
    r₁ = collect(lat_lon_to_cartesian(φ₁, λ₁))
    r₂ = collect(lat_lon_to_cartesian(φ₂, λ₂))

    n = cross(r₁, r₂)
    nrm = norm(n)
    if !(nrm > 0) || !isfinite(nrm)
        throw(ArgumentError("Great-circle normal is undefined for coincident or antipodal points."))
    end
    n̂ = n / nrm

    t₁ = cross(n̂, r₁); t₁ /= norm(t₁)
    t₂ = cross(n̂, r₂); t₂ /= norm(t₂)

    num = dot(n̂, cross(t₁, t₂))
    den = dot(t₁, t₂)

    return atan(num, den)
end

"""
    spherical_distance(a₁, a₂)
Compute the great-circle arc angle (in radians) between two points on the sphere, given their Cartesian coordinates `a₁`
and `a₂`. Both inputs are expected to be 3-vectors of same norm.

# Arguments
- `a₁`, `a₂`: 3-element Cartesian vectors on the sphere.

# Returns
- The spherical angle (in radians) between `a₁` and `a₂`.
"""
function spherical_distance(a₁, a₂)
    (sum(a₁.^2) ≈ sum(a₂.^2)) || error("a₁ and a₂ must have same norm")

    φ₁, λ₁ = cartesian_to_lat_lon(a₁)
    φ₂, λ₂ = cartesian_to_lat_lon(a₂)

    return haversine((λ₁, φ₁), (λ₂, φ₂), 1)
end

"""
    spherical_area_triangle(a::Number, b::Number, c::Number)

Returns the area of a spherical triangle on the unit sphere with sides `a`, `b`, and `c`.

The area of a spherical triangle on the unit sphere is ``E = A + B + C - π``, where ``A``, ``B``, and ``C`` are the
triangle's inner angles.

It has been known since the time of Euler and Lagrange that
``\\tan(E/2) = P / (1 + \\cos a + \\cos b + \\cos c)``, where
``P = (1 - \\cos²a - \\cos²b - \\cos²c + 2 \\cos a \\cos b \\cos c)^{1/2}``.

References
==========

* Euler, L. (1778) De mensura angulorum solidorum, Opera omnia, 26, 204-233 (Orig. in Acta adac. sc. Petrop. 1778)
* Lagrange,  J.-L. (1798) Solutions de quilquies problèmes relatifs au triangles sphéruques, Oeuvres, 7, 331-359.
"""
function spherical_area_triangle(a::Number, b::Number, c::Number)
    cosa = cos(a)
    cosb = cos(b)
    cosc = cos(c)

    tan½E = sqrt(1 - cosa^2 - cosb^2 - cosc^2 + 2cosa * cosb * cosc)
    tan½E /= 1 + cosa + cosb + cosc

    return 2atan(tan½E)
end

"""
    spherical_area_triangle(a₁, a₂, a₃)

Returns the area of a spherical triangle on the unit sphere with vertices given by the 3-vectors `a₁`, `a₂`, and `a₃`,
whose origin is the center of the sphere. The formula was first given by Eriksson (1990).

If we denote with ``A``, ``B``, and ``C`` the inner angles of the spherical triangle and with ``a``, ``b``, and ``c`` 
the sides of the triangle, then it has been known since Euler and Lagrange that
``\\tan(E/2) = P / (1 + \\cos a + \\cos b + \\cos c)``, where ``E = A + B + C - π`` is the triangle's excess and 
``P = (1 - \\cos²a - \\cos²b - \\cos²c + 2 \\cos a \\cos b \\cos c)^{1/2}``. 

On the unit sphere, ``E`` is precisely the area of the spherical triangle. Eriksson (1990) showed that ``P`` above is
the same as the volume defined by the vectors `a₁`, `a₂`, and `a₃`, that is ``P = |𝐚₁ ⋅ (𝐚₂ × 𝐚₃)|``.

References
==========

* Eriksson, F. (1990) On the measure of solid angles, Mathematics Magazine, 63 (3), 184-187, 
doi:10.1080/0025570X.1990.11977515
"""

function spherical_area_triangle(a₁, a₂, a₃)
    a₁, a₂, a₃ = collect(a₁), collect(a₂), collect(a₃)
    (sum(a₁.^2) ≈ 1 && sum(a₂.^2) ≈ 1 && sum(a₃.^2) ≈ 1) || error("a₁, a₂, a₃ must be unit vectors")

    tan½E = abs(dot(a₁, cross(a₂, a₃)))
    tan½E /= 1 + dot(a₁, a₂) + dot(a₂, a₃) + dot(a₁, a₃)

    return 2atan(tan½E)
end

"""
    spherical_area_quadrilateral(a₁, a₂, a₃, a₄)

Returns the area of a spherical quadrilateral on the unit sphere whose points are given by 3-vectors, `a₁`, `a₂`, `a₃`,
and `a₄`. The area of the quadrilateral is given as the sum of the areas of the two non-overlapping triangles. To avoid
having to pick the triangles appropriately ensuring they are not overlapping, we compute the area of the quadrilateral
as half the sum of the areas of all four potential triangles formed by `a₁`, `a₂`, `a₃`, and `a₄`.
"""
spherical_area_quadrilateral(a₁, a₂, a₃, a₄) =
    1/2 * (spherical_area_triangle(a₁, a₂, a₃) + spherical_area_triangle(a₁, a₂, a₄) +
           spherical_area_triangle(a₁, a₃, a₄) + spherical_area_triangle(a₂, a₃, a₄))

"""
    spherical_quadrilateral_vertices(X, Y, Z, i, j)

Returns the four Cartesian vertex vectors of the spherical grid cell whose corners are indexed by `(i, j)`, `(i+1, j)`,
`(i+1, j+1)`, and `(i, j+1)` in the arrays `X`, `Y`, and `Z`. Each of `X`, `Y`, and `Z` is a 2D array of size `(Nx, Ny)`
holding the Cartesian coordinates of grid vertices on the sphere, such that the point at `(i, j)` is
`(X[i, j], Y[i, j], Z[i, j])`.

# Arguments
- `X`, `Y`, `Z`: `(Nx, Ny)` arrays of Cartesian coordinates on the sphere.
- `i`, `j`: Indices of the lower-left corner of the cell (in array order).

# Returns
- `(a₁, a₂, a₃, a₄)`: The four 3-element Cartesian vertex vectors at `(i, j)`, `(i+1, j)`, `(i+1, j+1)`, and `(i, j+1)`.
"""
function spherical_quadrilateral_vertices(X, Y, Z, i, j)
    x₁ = X[i, j]
    y₁ = Y[i, j]
    z₁ = Z[i, j]
    a₁ = [x₁, y₁, z₁]
    x₂ = X[i+1, j]
    y₂ = Y[i+1, j]
    z₂ = Z[i+1, j]
    a₂ = [x₂, y₂, z₂]
    x₃ = X[i+1, j+1]
    y₃ = Y[i+1, j+1]
    z₃ = Z[i+1, j+1]
    a₃ = [x₃, y₃, z₃]
    x₄ = X[i, j+1]
    y₄ = Y[i, j+1]
    z₄ = Z[i, j+1]
    a₄ = [x₄, y₄, z₄]

    return a₁, a₂, a₃, a₄
end

"""
    compute_deviation_from_isotropy(X, Y, Z)

Compute a scalar measure of the deviation from isotropy for a spherical grid (e.g., a conformal cubed-sphere panel),
defined by the coordinate arrays `X`, `Y`, and `Z`. Each of `X`, `Y`, and `Z` is a 2D array of size `(Nx, Ny)` holding
the Cartesian coordinates of grid vertices on the sphere, such that the point at `(i, j)` is
`(X[i, j], Y[i, j], Z[i, j])`. The grid therefore contains `(Nx-1) × (Ny-1)` spherical quadrilateral cells.

For each quadrilateral cell in the grid, the function computes the arc lengths of its four edges on the unit sphere and
evaluates the sum of consecutive edge differences as a measure of anisotropy. The total deviation is then returned as 
the Euclidean norm over all cells.

# Arguments
- `X`, `Y`, `Z`: `(Nx, Ny)` arrays of Cartesian coordinates on the sphere.

# Returns
- A nonnegative scalar quantifying overall grid anisotropy (larger ⇒ more anisotropic).
"""
function compute_deviation_from_isotropy(X, Y, Z)
    Nx, Ny = size(X)
    deviation_from_isotropy = zeros(Nx-1, Ny-1)

    for j in 1:Ny-1, i in 1:Nx-1
        a₁, a₂, a₃, a₄ = spherical_quadrilateral_vertices(X, Y, Z, i, j)

        # Compute the arc lengths (distances) between the points a₁ and a₂, a₂ and a₃, a₃ and a₄, and a₄ and a₁ on the
        # unit sphere.
        d₁ = spherical_distance(a₁, a₂)
        d₂ = spherical_distance(a₂, a₃)
        d₃ = spherical_distance(a₃, a₄)
        d₄ = spherical_distance(a₄, a₁)

        # Compute the deviation from isotropy.
        deviation_from_isotropy[i, j] = abs(d₁ - d₂) + abs(d₂ - d₃) + abs(d₃ - d₄) + abs(d₄ - d₁)
    end

    return norm(deviation_from_isotropy)
end

"""
    compute_cell_areas(X, Y, Z)

Compute the spherical surface areas of all quadrilateral cells in a spherical grid (e.g., a conformal cubed sphere 
panel) defined by the coordinate arrays `X`, `Y`, and `Z`. Each of `X`, `Y`, and `Z` is a 2D array of size `(Nx, Ny)` 
holding the Cartesian coordinates of the grid vertices on the sphere, such that the point at `(i, j)` is
`(X[i, j], Y[i, j], Z[i, j])`. The grid therefore contains `(Nx-1) × (Ny-1)` spherical quadrilateral cells.

For each cell centered at `(i, j)` with vertices `(i, j)`, `(i+1, j)`, `(i+1, j+1)`, and `(i, j+1)`, the function
computes the cell area using `spherical_area_quadrilateral` and stores the results in a 2D array.

# Arguments
- `X`, `Y`, `Z`: `(Nx, Ny)` arrays of Cartesian coordinates on the sphere.

# Returns
- `cell_areas`: An `(Nx-1, Ny-1)` array of spherical quadrilateral cell areas.
"""
function compute_cell_areas(X, Y, Z)
    Nx, Ny = size(X)
    cell_areas = zeros(Nx-1, Ny-1)

    for j in 1:Ny-1, i in 1:Nx-1
        a₁, a₂, a₃, a₄ = spherical_quadrilateral_vertices(X, Y, Z, i, j)
        cell_areas[i, j] = spherical_area_quadrilateral(a₁, a₂, a₃, a₄)
    end

    return cell_areas
end
