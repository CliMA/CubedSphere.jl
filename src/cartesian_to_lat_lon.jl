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
    return atand(z, hypot(x, y))
end

function cartesian_to_longitude(x, y, z)
    return atand(y / x)
end
