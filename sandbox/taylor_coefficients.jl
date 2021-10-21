# Failed attempt to reproduce the Taylor coefficients in appendix B of Rančić et al. (1996).

function taylor_coefficients(A)
    A₁_expected = 1.47713062600964
    for _ in 1:10
        @printf("A₁ ≈ %.16f + %.16fim, ΔA₁ = %.16f + %.16fim\n", real(A[1]), imag(A[1]), real(A₁_expected - A[1]), imag(A₁_expected - A[1]))
        ϕ = range(-π/4, π/4, length=101)
        z = @. r * exp(im * ϕ)
        z′ = @. 1 - z
        Z = @. z′ .^ 4
        W = [sum(A[k] * Z[n]^k for k in 1:length(A)) for n in 1:length(Z)]
        g = [W[n] * r^4 * exp(im * θ) for (n, θ) in enumerate(range(-π, π, length=length(W)))]
        g̃ = fft(g)
        A = [g̃[k] / r^(4k) for k in 1:length(g̃)]
    end
end
