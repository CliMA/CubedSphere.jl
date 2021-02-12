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
    # xc=abs(x);
    # yc=abs(y);

    # kx=find(x<0);
    # ky=find(y<0);
    # kxy=find(yc>xc);

    # X=xc;
    # Y=yc;
    # xc=1-xc;
    # yc=1-yc;
    # xc(kxy)=1-Y(kxy);
    # yc(kxy)=1-X(kxy);

    # z=((xc+i*yc)/2).^4;
    # W=WofZ(z);

    X = xᶜ = abs(x)
    Y = yᶜ = abs(y)

    kxy = yᶜ > xᶜ

    xᶜ = 1 - xᶜ
    yᶜ = 1 - yᶜ

    kxy && (xᶜ = 1 - Y)
    kxy && (yᶜ = 1 - X)

    Z = ((xᶜ + im * yᶜ) / 2)^4
    W = W_Rancic(Z)

    # thrd=1/3;
    # i3=i^thrd;
    # ra=sqrt(3)-1;
    # cb=i-1;
    # cc=ra*cb/2;

    # W=i3*(W*i).^thrd;
    # W=(W-ra)./(cb+cc*W);
    # X=real(W);
    # Y=imag(W);
    # H=2./(1+X.^2+Y.^2);
    # X=X.*H;
    # Y=Y.*H;
    # Z=H-1;

    # T=X;
    # X(kxy)=Y(kxy);
    # Y(kxy)=T(kxy);

    # Y(ky)=-Y(ky);
    # X(kx)=-X(kx);

    # % Fix truncation for x=0 or y=0 (aja)
    # X(find(x==0))=0;
    # Y(find(y==0))=0;

    im³ = im^(1/3)
    ra = √3 - 1
    cb = -1 + im
    cc = ra * cb / 2

    W = im³ * (W * im)^(1/3)
    W = (W - ra) / (cb + cc * W)
    X, Y = reim(Z)

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

# using GLMakie

# x = range(-1, 1, length=51)
# y = range(-1, 1, length=51)
# X = zeros(length(x), length(y))
# Y = zeros(length(x), length(y))
# Z = zeros(length(x), length(y))

# for (i, x′) in enumerate(x), (j, y′) in enumerate(y)
#     X[i, j], Y[i, j], Z[i, j] = conformal_cubed_sphere_mapping(x′, y′)
# end

# # wireframe(X, Y, Z)
# # surface(X, Y, Z)
# scatter(X[:], Y[:], Z[:])
