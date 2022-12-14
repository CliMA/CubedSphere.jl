using FFTW, SpecialFunctions, TaylorSeries, ProgressBars

function find_angles(φ)
    φ⁻ = -φ
    φ⁺ = +φ

    w⁻ = cis(φ⁻)
    w⁺ = cis(φ⁺)

    w′⁻ = (1 - w⁻) / (1 + w⁻/2)
    w′⁺ = (1 - w⁺) / (1 + w⁺/2)

    φ′⁻ = angle(w′⁻)
    φ′⁺ = angle(w′⁺)

    return φ′⁻, φ′⁺
end

# cbrt goes from W to w
# cbrt′ goes from W′ to w′
function Base.cbrt(z::Complex)
    r = abs(z)
    φ = angle(z) # ∈ [-π, +π]
    θ = φ / 3
    return cbrt(r) * cis(θ)
end

function cbrt′(z::Complex)
    φ′⁻, φ′⁺ = find_angles(π/3)
    r = abs(z)
    φ = angle(z) # ∈ [-π, +π]

    θ = φ / 3

    if 0 < θ ≤ φ′⁻
        θ -= 2π/3
    elseif φ′⁺ ≤ θ < 0
        θ += 2π/3
    end
    return r^(1/3) * cis(θ)
end

"""
    find_N(r; decimals=15)

Return the required number of points we need to consider around the circle of
radius `r` to compute the conformal map series coefficients up to `decimals`
points. The number of points is computed based on the estimate of eq. (B9) in
the paper by Rančić et al. (1996).

Note: the returned number of points is a power of 2 so that the FFTs are efficient.
"""
function find_N(r; decimals=15)
    A₁ = 1.4771 # an approximation of the first coefficient
    C = sqrt(3) * gamma(1/3) * A₁^(1/3) / ((256)^(1/3) * π)

    N = 2
    while N + 7/12 * log10(N) / (-log10(r)) - (decimals + log10(A₁ / C)) / (-4 * log10(r)) < 0
        N *= 2
    end

    return N
end

function update_coefficients!(A, r, Nφ)
    Ncoefficients = length(A)

    Lφ = π/2
    dφ = Lφ / Nφ
    φ = range(-Lφ/2 + dφ/2, stop=Lφ/2 - dφ/2, length=Nφ)

    z = @. r * exp(im * φ)

    W̃′ = 0z
    for k = 1:Ncoefficients
        @. W̃′ += A[k] * (1 - z)^(4k)
    end

    w̃′ = @. cbrt′(W̃′)
    w̃  = @. (1 - w̃′) / (1 + w̃′/2)
    W̃  = @. w̃^3

    k = collect(fftfreq(Nφ, Nφ))
    g̃ = fft(W̃) ./ (Nφ * exp.(im * k * 4φ[1]))
    g̃ = g̃[2:Ncoefficients+1] # exclude coefficient for 0th power

    A = [real(g̃[k] / r^(4k)) for k in 1:Ncoefficients]

    return nothing
end

"""
    find_taylor_coefficients(r = 1 - 1e-7; maximum_coefficients=256, Niterations=30)

Return the Taylor coefficients for the conformal map ``Z \to W`` and its inverse,
where ``Z = z^4`` and ``W = w^3``. In particular, it returns the coefficients
of the Taylor series
```math
W(Z) = \\sum_{k=1}^\\infty A_k Z^k
```

and its inverse

```math
Z(W) = \\sum_{k=1}^\\infty A_k Z^k
```

The algorithm to obtain the coefficients follows the procedure described by
Rančić et al. (1996) in their Appendix B.

References
==========

- Rančić et al., (1996). A global shallow-water model using an expanded spherical cube - Gnomonic versus conformal
  coordinates, _Quarterly Journal of the Royal Meteorological Society_.
"""
function find_taylor_coefficients(r = 1 - 1e-7; maximum_coefficients=256, Niterations=30)
    r ≥ 1 && error("r has to be less than 1")

    Nφ = find_N(r; decimals=15)
    Ncoefficients = Int(Nφ/2) - 2 > maximum_coefficients ? maximum_coefficients : Int(Nφ/2) - 2

    @info "Computing the first $Ncoefficients coefficients in the Taylor series"
    @info "using $Nφ function evaluations on a circle with radius $r."

    # initialize coefficients
    A_coefficients = rand(Ncoefficients)

    for _ in ProgressBar(1:Niterations)
        global A
        update_coefficients!(A_coefficients, r, Nφ)
    end

    # convert to Taylor series; add coefficient for 0-th power
    A_series = Taylor1([0; A_coefficients])

    B_series = inverse(A_series) # This is the inverse Taylor series

    B_series.coeffs[1] !== 0.0 && error("coefficient that corresponds to W^0 is non-zero; something went wrong")
    B_coefficients = B_series.coeffs[2:end]

    return A_coefficients, B_coefficients
end
