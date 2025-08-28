"""
    spherical_distance(a₁, a₂)

Compute the great-circle arc length (angle in radians) between two points on the **unit sphere**, given their Cartesian
coordinates `a₁` and `a₂`. Both inputs are expected to be 3-vectors of unit length; the dot product is clamped to
`[-1, 1]` to guard against floating-point roundoff before applying `acos`.

# Arguments
- `a₁`, `a₂`: 3-element Cartesian vectors on the unit sphere (‖`a`‖ = 1).

# Returns
- The arc length (in radians) between `a₁` and `a₂`.
"""
function spherical_distance(a₁::AbstractVector, a₂::AbstractVector)
    (sum(a₁.^2) ≈ 1 && sum(a₂.^2) ≈ 1) || error("a₁ and a₂ must be unit vectors")

    # Compute the dot product and calculate the arccosine to find the angle.
    cosθ = dot(a₁, a₂)

    # Ensure the result is within the domain of acos due to potential floating-point errors.
    cosθ = clamp(cosθ, -1, 1)

    # Return the arc length, which is the angle between the two points.
    return acos(cosθ)
end

"""
    spherical_area_triangle(a::Number, b::Number, c::Number)

Return the area of a spherical triangle on the unit sphere with sides `a`, `b`, and `c`.

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
    spherical_area_triangle(a::AbstractVector, b::AbstractVector, c::AbstractVector)

Return the area of a spherical triangle on the unit sphere with vertices given by the 3-vectors `a`, `b`, and `c`
whose origin is the the center of the sphere. The formula was first given by Eriksson (1990).

If we denote with ``A``, ``B``, and ``C`` the inner angles of the spherical triangle and with ``a``, ``b``, and ``c`` 
the side of the triangle, then it has been known since Euler and Lagrange that 
``\\tan(E/2) = P / (1 + \\cos a + \\cos b + \\cos c)``, where ``E = A + B + C - π`` is the triangle's excess and 
``P = (1 - \\cos²a - \\cos²b - \\cos²c + 2 \\cos a \\cos b \\cos c)^{1/2}``. 

On the unit sphere, ``E`` is precisely the area of the spherical triangle. Erikkson (1990) showed that ``P`` above is 
the same as the volume defined by the vectors `a`, `b`, and `c`, that is ``P = |𝐚 \\cdot (𝐛 \\times 𝐜)|``.

References
==========

* Eriksson, F. (1990) On the measure of solid angles, Mathematics Magazine, 63 (3), 184-187, 
doi:10.1080/0025570X.1990.11977515
"""
function spherical_area_triangle(a₁::AbstractVector, a₂::AbstractVector, a₃::AbstractVector)
    (sum(a₁.^2) ≈ 1 && sum(a₂.^2) ≈ 1 && sum(a₃.^2) ≈ 1) || error("a₁, a₂, a₃ must be unit vectors")

    tan½E = abs(dot(a₁, cross(a₂, a₃)))
    tan½E /= 1 + dot(a₁, a₂) + dot(a₂, a₃) + dot(a₁, a₃)

    return 2atan(tan½E)
end

"""
    spherical_area_quadrilateral(a₁, a₂, a₃, a₄)

Return the area of a spherical quadrilateral on the unit sphere whose points are given by 3-vectors, `a`, `b`, `c`, and
`d`. The area of the quadrilateral is given as the sum of the ares of the two non-overlapping triangles. To avoid having
to pick the triangles appropriately ensuring they are not overlapping, we compute the area of the quadrilateral as the
half the sum of the areas of all four potential triangles formed by `a₁`, `a₂`, `a₃`, and `a₄`.
"""
spherical_area_quadrilateral(a::AbstractVector, b::AbstractVector, c::AbstractVector, d::AbstractVector) =
    1/2 * (spherical_area_triangle(a, b, c) + spherical_area_triangle(a, b, d) +
           spherical_area_triangle(a, c, d) + spherical_area_triangle(b, c, d))

"""
    spherical_quadrilateral_vertices(X, Y, Z, i, j)

Return the four Cartesian vertex vectors of the spherical grid cell whose corners are indexed by `(i, j)`, `(i+1, j)`,
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
