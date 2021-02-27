"""
    conformal_cubed_sphere_mapping(x, y)

Conformal mapping of a cube onto a sphere. Maps `(x, y)` on the north-pole face of a cube
to (X, Y, Z) coordinates in physical space. The face is oriented normal to Z-axis with
X and Y increasing with x and y.

The input coordinates must lie within the range -1 <= x <= 1,  -1 <= y <= 1.

This numerical conformal mapping is described by Rančić et al. (1996).

This is a Julia translation of MATLAB code from MITgcm [1] that is based on
Fortran 77 code from Jim Purser & Misha Rančić.

[1] http://wwwcvs.mitgcm.org/viewvc/MITgcm/MITgcm_contrib/high_res_cube/matlab-grid-generator/map_xy2xyz.m?view=markup
"""
function conformal_cubed_sphere_mapping(x, y)
    X = xᶜ = abs(x)
    Y = yᶜ = abs(y)

    kxy = yᶜ > xᶜ

    xᶜ = 1 - xᶜ
    yᶜ = 1 - yᶜ

    kxy && (xᶜ = 1 - Y)
    kxy && (yᶜ = 1 - X)

    Z = ((xᶜ + im * yᶜ) / 2)^4
    W = W_Rancic(Z)

    im³ = im^(1/3)
    ra = √3 - 1
    cb = -1 + im
    cc = ra * cb / 2

    W = im³ * (W * im)^(1/3)
    W = (W - ra) / (cb + cc * W)
    X, Y = reim(W)

    H = 2 / (1 + X^2 + Y^2)
    X = X * H
    Y = Y * H
    Z = H - 1

    if kxy
        X, Y = Y, X
    end

    y < 0 && (Y = -Y)
    x < 0 && (X = -X)

    # Fix truncation for x = 0 or y = 0.
    x == 0 && (X = 0)
    y == 0 && (Y = 0)

    return X, Y, Z
end

W_Rancic(Z) = sum(A_Rancic[k] * Z^(k-1) for k in 1:length(A_Rancic))
