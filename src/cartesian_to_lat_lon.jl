"""
    cartesian_to_lat_lon(x, y, z)

Convert 3D coordinates `(x, y, z)` on the unit sphere to latitude-longitude `(lat, lon)`. Assumes `lat` is positive
with `z`, equatorial plane falls at ``z = 0``  and `lon` is measured anti-clockwise (eastward) from ``x``-axis (``y = 0``)
about the ``z``-axis.

This is a Julia translation of [MATLAB code from MITgcm](http://wwwcvs.mitgcm.org/viewvc/MITgcm/MITgcm_contrib/high_res_cube/matlab-grid-generator/map_xyz2lonlat.m?view=markup).

Examples
========

```@jlddoctest
julia> using CubedSphere

julia> cartesian_to_lat_lon(0, 0, 1)
(90.0, 0.0)

julia> cartesian_to_lat_lon(√2/4, -√2/4, √3/2)
(59.99999999999999, -45.0)
```
"""
cartesian_to_lat_lon(x, y, z) = cartesian_to_latitude(x, y, z), cartesian_to_longitude(x, y, z)

cartesian_to_latitude(x, y, z) = atand(z, hypot(x, y))

cartesian_to_longitude(x, y, z) = atand(y, x)
