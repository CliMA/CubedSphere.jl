"""
    cartesian_to_lat_lon(x, y, z)

Convert 3D cartesian coordinates `(x, y, z)` on the sphere to latitude-longitude. Returns a tuple
`(latitude, longitude)` in degrees.

The equatorial plane falls at ``z = 0``, latitude is the angle measured from the equatorial plane,
and longitude is measured anti-clockwise (eastward) from ``x``-axis (``y = 0``) about the ``z``-axis.

Examples
========

Find latitude-longitude of the North Pole

```jldoctest 1
julia> using CubedSphere

julia> x, y, z = (0, 0, 6.4e6); # cartesian coordinates of North Pole [in meters]

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

cartesian_to_latitude(x, y, z) = atand(z, hypot(x, y))

cartesian_to_longitude(x, y, z) = atand(y, x)
