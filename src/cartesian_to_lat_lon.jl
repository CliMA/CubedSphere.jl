"""
    cartesian_to_lat_lon(x, y, z)

Convert 3D coordinates (x, y, z) from the unit sphere to (lat, lon). Assumes "lat" is positive
with "z", equatorial plane falls at z=0  and "lon" is measured anti-clockwise (eastward)
from x-axis (y=0) about z-axis.

This is a Julia translation of MATLAB code from MITgcm [1].

[1]: http://wwwcvs.mitgcm.org/viewvc/MITgcm/MITgcm_contrib/high_res_cube/matlab-grid-generator/map_xyz2lonlat.m?view=markup
"""
cartesian_to_lat_lon(x, y, z) = cartesian_to_latitude(x, y, z), cartesian_to_longitude(x, y, z)

function cartesian_to_latitude(x, y, z)
    R_eq = √(x^2 + y^2)
    R_eq == 0 && return atand(z * Inf) # You're at one of the poles!
    return atand(z / R_eq)
end

function cartesian_to_longitude(x, y, z)
    x == 0 && return atand(y * Inf) # You're at the prime-meridian or anti-meridian!
    λ = atand(y / x)
    x < 0 && y >= 0 && return 180 + λ
    x <= 0 && y < 0 && return λ - 180
    return λ
end
