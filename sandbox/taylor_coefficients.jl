using GLMakie, Printf, FFTW, ProgressBars, SpecialFunctions, CubedSphere

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

Note: the returned number of points is a power of 2 so that
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

function update_coefficients(A, r, Nφ)
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
    g̃ = g̃[2:Ncoefficients+1] #exclude coefficient for z⁰

    A_update = [real(g̃[k] / r^(4k)) for k in 1:Ncoefficients]

    return A_update
end

function plot_transformation(A, r, Nφ; Lφ=π/2)
    dφ = Lφ / Nφ
    φ = range(-Lφ/2 + dφ/2, stop=Lφ/2 - dφ/2, length=Nφ)

    z = @. r * exp(im * φ)
    Z = @. z^4

    W  = zeros(eltype(z), size(z))
    W̃′ = zeros(eltype(z), size(z))
    for k = 1:length(A)
        @. W  += A[k] * z^(4k)
        @. W̃′ += A[k] * (1 - z)^(4k)
    end

    w = @. cbrt(W)
    w′ = @. (1 - w) / (1 + w/2)
    W′ = @. w′^3

    w̃′ = @. cbrt′(W̃′, r)
    w̃  = @. (1 - w̃′) / (1 + w̃′/2)
    W̃  = @. w̃^3

    fig = Figure(resolution=(1200, 1200), fontsize=30)

    axz = Axis(fig[1, 1], title="z")
    axZ = Axis(fig[1, 2], title="Z")
    axw = Axis(fig[1, 1], title="w")
    axW = Axis(fig[1, 2], title="W")
    axwp = Axis(fig[2, 1], title="w′")
    axWp = Axis(fig[2, 2], title="W′")

    lim = 2
    for ax in (axw, axW)
        xlims!(ax, -lim, lim)
        ylims!(ax, -lim, lim)
    end

    lim = 1
    for ax in (axwp, axWp)
        xlims!(ax, -lim, lim)
        ylims!(ax, -lim, lim)
    end

    lim = 1.5
    for ax in (axz, axZ)
        xlims!(ax, -lim, lim)
        ylims!(ax, -lim, lim)
    end

    scatter!(axz, real.(z), imag.(z), linewidth=4)
    scatter!(axZ, real.(Z), imag.(Z))

    scatter!(axw, real.(w), imag.(w), color=(:black, 0.5), linewidth=8, label=L"w")
    scatter!(axw, real.(w̃), imag.(w̃), color=(:orange, 0.8), linewidth=4, label=L"w̃")

    scatter!(axW, real.(W), imag.(W), color=(:black, 0.5), linewidth=8, label=L"W")
    scatter!(axW, real.(W̃), imag.(W̃), color=(:orange, 0.8), linewidth=4, label=L"W̃")

    scatter!(axwp, real.(w′), imag.(w′), color=(:black, 0.5), linewidth=8, label=L"w′")
    scatter!(axwp, real.(w̃′), imag.(w̃′), color=(:orange, 0.8), linewidth=4, label=L"w̃′")

    scatter!(axWp, real.(W′), imag.(W′), color=(:black, 0.5), linewidth=8, label=L"W′")
    scatter!(axWp, real.(W̃′), imag.(W̃′), color=(:orange, 0.8), linewidth=4, label=L"W̃′")

    for ax in [axw, axW, axwp, axWp]
        axislegend(ax)
    end

    save("consistency.png", fig)

    display(fig)
end

r = 1 - 1e-7
Nφ = find_N(r; decimals=15)

maximum_coefficients = 512
Ncoefficients = Int(Nφ/2) - 2 > maximum_coefficients ? maximum_coefficients : Int(Nφ/2) - 2

@info "Attempting to compute the first $Ncoefficients coefficients in the Taylor series."
@info "To do so, we make function evaluations on a circle with radius $r."

# initialize coefficients
A = rand(Ncoefficients)

Niterations = 30
for _ in ProgressBar(1:Niterations)
    global A
    A = update_coefficients(A, r, Nφ)
end

@info "After $Niterations iterations we have:"

for (k, Aₖ) in enumerate(A[1:30])
    @printf("k = %2i, A ≈ %+.14f, A_Rancic = %+.14f, |A - A_Rancic| = %.2e \n", k, real(Aₖ), CubedSphere.A_Rancic[k+1], abs(CubedSphere.A_Rancic[k+1] - real(Aₖ)))
end

plot_transformation(A, r, Nφ)
