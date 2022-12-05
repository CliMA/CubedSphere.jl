import Elliptic.Jacobi: sn, cn, dn

"""
    sn(z::Complex, m::Real)

Compute the Jacobi elliptic function `sn(z | m)` following Abramowitz & Stegun (1964), Eq. 16.21.2.
"""
@inline function sn(z::Complex, m::Real)
    s  = sn(real(z), m)
    s₁ = sn(imag(z), 1-m)
    c  = cn(real(z), m)
    c₁ = cn(imag(z), 1-m)
    d  = dn(real(z), m)
    d₁ = dn(imag(z), 1-m)
    return (s * d₁ + im * c * d * s₁ * c₁) / (c₁^2 + m * s^2 * s₁^2)
end

"""
    cn(z::Complex, m::Real)

Compute the Jacobi elliptic function `cn(z | m)` following Abramowitz & Stegun (1964), Eq. 16.21.3.
"""
@inline function cn(z::Complex, m::Real)
    s  = sn(real(z), m)
    s₁ = sn(imag(z), 1-m)
    c  = cn(real(z), m)
    c₁ = cn(imag(z), 1-m)
    d  = dn(real(z), m)
    d₁ = dn(imag(z), 1-m)
    return (c * c₁ - im * s * d * s₁ * d₁) / (c₁^2 + m * s^2 * s₁^2)
end

"""
    dn(z::Complex, m::Real)

Compute the Jacobi elliptic function `dn(z | m)` following Abramowitz & Stegun (1964), Eq. 16.21.4.
"""
@inline function dn(z::Complex, m::Real)
    s  = sn(real(z), m)
    s₁ = sn(imag(z), 1-m)
    c  = cn(real(z), m)
    c₁ = cn(imag(z), 1-m)
    d  = dn(real(z), m)
    d₁ = dn(imag(z), 1-m)
    return (d * c₁ * d₁ - im * m * s * c * s₁) / (c₁^2 + m * s^2 * s₁^2)
end
