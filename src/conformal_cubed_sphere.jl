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

"""
   lamrofnoc_cubed_sphere_mapping(X, Y, Z)

Inverse mapping for conformal cube sphere for quadrant of north-pole face in which X and Y are both
positive. All other mappings to other cube face coordinates can be recovered from rotations of this map.
There a 3 other quadrants for the north-pole face and five other faces for a total of twenty-four quadrants.
Because of symmetry only the reverse for a single quadrant is needed. Because of branch cuts and the complex
transform the inverse mappings are multi-valued in general, using a single quadrant case allows a simple
set of rules to be applied.

"""
function lamrofnoc_cubed_sphere_mapping(X, Y, Z)

        H  = Z + 1
        Xˢ = X/H
        Yˢ = Y/H
        ω  = Xˢ + im*Yˢ

        ra = √3 - 1
        cb = -1 + im
        cc = ra * cb / 2
        ω⁰ = (ω*cb + ra)/(1-ω*cc)
        W⁰ = im*ω⁰^3*im
        Z  = Z_Rancic(W⁰)
        z  = (Z^0.25)*2
        x, y = reim(z)

        kxy = abs(y) > abs(x)
        xx = x
        yy = y
        !kxy && ( x = 1 - abs(yy) )
        !kxy && ( y = 1 - abs(xx) )

        xf = x
        yf = y
        ( X < Y ) && ( xf = y  )
        ( X < Y ) && ( yf = x  )
        x = xf
        y = yf

        return x, y
end

W_Rancic(Z) = sum(A_Rancic[k] * Z^(k-1) for k in 1:length(A_Rancic))
Z_Rancic(W) = sum(B_Rancic[k] * W^(k-1) for k in 1:length(B_Rancic))
