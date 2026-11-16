using InverseProblems
using LinearAlgebra
using Random
using SparseArrays
using Test

# Simply supported beam, modal synthesis of the transfer functions
function beam_frf(ω; nm = 30, nr = 25, nmodes = 30, ξ = 0.01)
    xm = range(0, 1, length = nm + 2)[2:end-1]
    xr = range(0, 1, length = nr + 2)[2:end-1]
    ωn = [(n * π)^2 * 50 for n in 1:nmodes]
    return [sum(sin(n * π * xm[i]) * sin(n * π * xr[j]) / (ωn[n]^2 - ω^2 + 2im * ξ * ωn[n] * ω)
                for n in 1:nmodes) for i in 1:nm, j in 1:nr]
end

const ωn5, ωn6 = (5π)^2 * 50, (6π)^2 * 50

@testset "InverseProblems.jl" begin
    @testset "Problem construction" begin
        H = randn(10, 6)
        @test size(InverseProblem(H, randn(10)).y) == (10,)
        @test_throws DimensionMismatch InverseProblem(H, randn(9))
        # space–time: H acts on vec(U), U is 3 × 2
        p = InverseProblem(randn(10, 6), randn(5, 2))
        @test InverseProblems.unknown_size(p) == (3, 2)
        # space–frequency, vector of blocks or 3D array
        Hs = [randn(4, 3) for _ in 1:5]
        @test InverseProblems.unknown_size(InverseProblem(Hs, randn(4, 5))) == (3, 5)
        @test InverseProblems.unknown_size(InverseProblem(randn(4, 3, 5), randn(4, 5))) == (3, 5)
    end

    @testset "Regularization terms" begin
        @test_throws ArgumentError Lq(0.0)
        @test_throws ArgumentError Lq(3.0)
        @test_throws ArgumentError Lpq(2, 2.5)
        U = [1.0 -2.0; 0.0 3.0]
        @test InverseProblems.penalty(Lq(1.0), U) ≈ 6
        @test InverseProblems.penalty(Lpq(2, 1), U) ≈ sqrt(5) + 3
        @test InverseProblems.penalty(Lpq(2, 2), U) ≈ InverseProblems.penalty(Lq(2.0), U)
        G = Groups(1:1 => Lq(1.0), 2:2 => Lq(2.0))
        R = InverseProblems.Regularization(G, 2)
        @test R.L isa InverseProblems.IdentitySmoothing
        @test InverseProblems.penalty(R, U) ≈ 3 + 9
        # one smoothing matrix per group, assembled into a global operator
        D = [1.0 -1.0]
        R = InverseProblems.Regularization(Groups([1, 2] => Lq(2.0; L = D), [3] => Lq(1.0)), 3)
        @test Matrix(R.L.L) == [1 -1 0; 0 0 1]
        @test R.rows == [[1], [2]]
        @test_throws DimensionMismatch InverseProblems.groups(Groups([1] => Lq(2.0; L = D), [2] => Lq(1.0)), 2)
        @test_throws ArgumentError InverseProblems.groups(Groups(1:2 => Lq(1.0), 2:3 => Lq(2.0)), 3)
        @test_throws ArgumentError InverseProblems.groups(Groups(1:2 => Lq(1.0)), 3)
    end

    @testset "Formulations" begin
        @test_throws ArgumentError Additive(Lq(2.0); λ = -1.0)
        @test_throws ArgumentError Additive(Lq(2.0); λ = Heuristic())
        @test InverseProblems.default_algorithm(Multiplicative(Lq(1.0))) isa IR
        @test InverseProblems.default_algorithm(Bayesian(Lq(1.0))) isa MAP
    end

    @testset "Weighted solves" begin
        Random.seed!(0)
        IP = InverseProblems
        # every regime × smoothing operator × field, against the explicit formulas
        for T in (Float64, ComplexF64), (nm, nr) in ((12, 8), (8, 8), (5, 8)), ns in (1, 3)
            Lsq = diagm(nr, nr, 0 => fill(2.0, nr), 1 => -ones(nr - 1))       # square invertible
            Lrect = diagm(nr, nr, 0 => -ones(nr), 1 => ones(nr - 1))[1:nr-1, :]   # rectangular
            for (Lop, Lm) in ((IP.IdentitySmoothing(), Matrix(1.0I, nr, nr)),
                              (IP.MatrixSmoothing(sparse(Lsq)), Lsq),
                              (IP.MatrixSmoothing(sparse(Lrect)), Lrect))
                Lg = kron(I(ns), Lm)
                # dense operator acting on vec(U) (space–time layout when ns > 1)
                H, y = randn(T, nm * ns, nr * ns), randn(T, nm, ns)
                c = IP.linear_cache(InverseProblem(H, ns == 1 ? vec(y) : y))
                @test (c isa IP.Underdetermined) == (nm < nr)
                D = rand(size(Lm, 1), ns) .+ 0.1
                A = H'H + Lg' * Diagonal(vec(D)) * Lg
                @test vec(IP.weighted_solve(c, Lop, D)) ≈ A \ (H' * vec(y))
                @test IP.posterior_covariance(c, Lop, D) ≈ inv(A)
                # block-diagonal operator
                Hs = [randn(T, nm, nr) for _ in 1:ns]
                cb = IP.linear_cache(InverseProblem(Hs, y))
                U = IP.weighted_solve(cb, Lop, D)
                Σ = IP.posterior_covariance(cb, Lop, D)
                for j in 1:ns
                    Aj = Hs[j]'Hs[j] + Lm' * Diagonal(D[:, j]) * Lm
                    @test U[:, j] ≈ Aj \ (Hs[j]' * y[:, j])
                    @test Σ[j] ≈ inv(Aj)
                end
            end
        end
        @test IP.MatrixSmoothing(sparse([1.0 1; 1 1])) isa IP.SingularSmoothing
        @test IP.MatrixSmoothing(sparse([1.0 1; 0 1])) isa IP.InvertibleSmoothing
        # Tikhonov with a fixed parameter: closed form
        H, y = randn(15, 10), randn(15)
        s = solve(InverseProblem(H, y), Additive(Lq(2.0); λ = 0.3))
        @test s.u ≈ (H'H + 0.3I) \ (H'y)
        @test s.converged
        # ... with a smoothing operator (first-order differences, non-square)
        Dx = diagm(10, 10, 0 => -ones(10), 1 => ones(9))[1:9, :]
        s = solve(InverseProblem(H, y), Additive(Lq(2.0; L = Dx); λ = 0.3))
        @test s.u ≈ (H'H + 0.3Dx'Dx) \ (H'y)
    end

    @testset "Parameter choice" begin
        Random.seed!(3)
        H = randn(20, 12)
        y = H * randn(12) .+ 0.1randn(20)
        c = InverseProblems.linear_cache(InverseProblem(H, y))
        W = rand(12, 1) .+ 0.5
        Lid = InverseProblems.IdentitySmoothing()
        Lmat = InverseProblems.MatrixSmoothing(sparse(1.0I, 12, 12))
        # the SVD and generalized SVD paths describe the same problem
        for rule in (LCurve(), GCV(), BayesEstimator())
            λ1 = InverseProblems.select_parameter(rule, c, Lid, W)
            λ2 = InverseProblems.select_parameter(rule, c, Lmat, W)
            @test λ1 > 0
            @test isapprox(λ1, λ2; rtol = 1e-6)
        end
        # 1D refinement: Brent's method, and refined grid search better than the grid alone
        x, fx = InverseProblems.brent(x -> (x - 0.3)^2 + 1, -2, 3)
        @test x ≈ 0.3 atol = 1e-7
        @test fx ≈ 1
        @test InverseProblems.grid_minimize(x -> (x - 0.3)^2 + 1, range(-2, 3; length = 11)) ≈ 0.3 atol = 1e-7
        sp = InverseProblems.spectrum(c, Lid, W)
        t = sort(log.(InverseProblems.lambda_grid(sp.γ)))
        for rule in (LCurve(), GCV(), BayesEstimator())
            J = InverseProblems.criterion(rule, sp, length(y))
            λ = InverseProblems.select_parameter(rule, sp, length(y))
            @test J(log(λ) / 2) <= minimum(J.(t))
        end
        s = solve(InverseProblem(H, y), Additive(Lq(2.0); λ = BayesEstimator()))
        @test s.λ > 0 && s.converged
        # non-square L: the null space of L is left unregularized
        Dx = diagm(12, 12, 0 => -ones(12), 1 => ones(11))[1:11, :]
        sp = InverseProblems.spectrum(c, InverseProblems.MatrixSmoothing(sparse(Dx)), ones(11, 1))
        @test sp.nfree == 1 && length(sp.γ) == 11
    end

    @testset "Point force reconstruction (single frequency)" begin
        Random.seed!(1)
        H = beam_frf((ωn5 + ωn6) / 2)
        u = zeros(ComplexF64, size(H, 2))
        u[17] = 1
        y0 = H * u
        prob = InverseProblem(H, y0 .+ 0.05norm(y0) / sqrt(length(y0)) .* randn(ComplexF64, length(y0)))
        sparse = Groups(1:12 => Lq(1.0), 13:25 => Lq(0.5))
        cases = [
            (Additive(Lq(1.0); λ = GCV()), IR()),
            (Multiplicative(Lq(0.5)), IR()),
            (Additive(Lq(1.0); λ = BayesEstimator()), IR()),
            (Multiplicative(sparse; multiparameter = true), IR()),
            (Bayesian(Lq(1.0)), MAP()),
        ]
        for (reg, alg) in cases
            s = solve(prob, reg, alg)
            @test argmax(abs.(s.u)) == 17
            @test norm(s.u - u) / norm(u) < 0.1
        end
        s = solve(prob, Bayesian(Lq(1.0)))
        @test all(0.05 .<= s.hyper.q .<= 2)
        @test s.hyper.τn > 0
    end

    @testset "Bayesian, real data" begin
        Random.seed!(4)
        H = real(beam_frf((ωn5 + ωn6) / 2))
        u = zeros(25)
        u[9] = 1
        y0 = H * u
        σ = 0.02norm(y0) / sqrt(length(y0))
        s = solve(InverseProblem(H, y0 .+ σ .* randn(length(y0))), Bayesian(Lq(1.0)))
        @test eltype(s.u) <: Real
        @test argmax(abs.(s.u)) == 9
        @test norm(s.u - u) / norm(u) < 0.1
        # real model: τn estimates 1/σ²
        @test 0.1 < s.hyper.τn * σ^2 < 3
        # λ_g = τs_g / τn (MAP ½‖r‖² + (λ/q)‖u‖_q^q in the real case)
        @test s.λ ≈ s.hyper.τs ./ s.hyper.τn
    end

    @testset "Space–frequency mixed norm" begin
        Random.seed!(2)
        ωs = range(300, 15000, length = 20)
        Hs = [beam_frf(ω) for ω in ωs]
        U = zeros(ComplexF64, 25, length(ωs))
        U[17, :] .= 1
        Y = reduce(hcat, [Hs[j] * U[:, j] for j in eachindex(ωs)])
        prob = InverseProblem(Hs, Y)
        s = solve(prob, Multiplicative(Lpq(2, 0.5)))
        e = vec(sum(abs2, s.u; dims = 2))
        @test argmax(e) == 17
        @test e[17] / sum(e) > 0.95
        @test size(s.u) == size(U)
        # posterior covariance: one nr × nr matrix per frequency
        s = solve(prob, Additive(Lq(2.0); λ = 1e-12); covariance = true)
        @test length(s.covariance) == length(ωs)
        @test size(posterior_std(s)) == size(U)
        @test size(posterior_samples(s, 10)) == (size(U)..., 10)
        lo, up = credible_interval(s; nsamples = 500)
        @test size(lo) == size(U) && all(real(lo) .<= real(s.u) .<= real(up))
    end

    @testset "Noise covariance" begin
        Random.seed!(6)
        H, y = randn(15, 10), randn(15)
        prob = InverseProblem(H, y)
        # σ²I: equivalent to a rescaled regularization parameter
        s = solve(prob, Additive(Lq(2.0); λ = 0.3); R = 4.0I)
        @test s.u ≈ (H'H + 1.2I) \ (H'y)
        # diagonal R given as a vector, a Diagonal or a full diagonal matrix
        v = rand(15) .+ 0.5
        ref = (H' * (H ./ v) + 0.3I) \ (H' * (y ./ v))
        for R in (v, Diagonal(v), Matrix(Diagonal(v)))
            @test solve(prob, Additive(Lq(2.0); λ = 0.3); R).u ≈ ref
        end
        @test_throws ArgumentError solve(prob, Additive(Lq(2.0); λ = 0.3); R = ones(15, 15))
        @test_throws DimensionMismatch solve(prob, Additive(Lq(2.0); λ = 0.3); R = ones(14))
        # posterior covariance of the additive formulation
        s = solve(prob, Additive(Lq(2.0); λ = 0.3); R = v, covariance = true)
        @test s.covariance ≈ inv(H' * (H ./ v) + 0.3I)
        @test posterior_std(s) ≈ sqrt.(diag(inv(H' * (H ./ v) + 0.3I)))
        @test_throws ArgumentError posterior_std(solve(prob, Additive(Lq(2.0); λ = 0.3)))
    end

    @testset "Generalized Gaussian priors (real / complex)" begin
        IP = InverseProblems
        # real prior ∝ exp(-(τs/q)|v|^q), complex prior ∝ exp(-τs|v|^q): the alternating
        # estimation of (q, τs) recovers the parameters of the sampled distribution
        Random.seed!(8)
        n, τ = 20_000, 4.0
        cases = [(IP.RealField(), 2.0, randn(n) / sqrt(τ)),                          # N(0, 1/τ)
                 (IP.RealField(), 1.0, (2 .* rand(Bool, n) .- 1) .* randexp(n) ./ τ),  # Laplace
                 (IP.ComplexField(), 2.0, randn(ComplexF64, n) / sqrt(τ))]            # CN(0, 1/τ)
        for (F, qtrue, v) in cases
            q, τs = 1.0, 1.0
            for _ in 1:50
                q = IP.optimal_shape(F, v, τs, (0.05, 2.0))
                τs = IP.signal_precision(F, v, q)
            end
            @test q ≈ qtrue rtol = 0.05
            @test τs ≈ τ rtol = 0.1
        end
    end

    @testset "Element-wise Bayesian (RVR)" begin
        # reference: direct port of the original RVR implementation
        function rvr_ref(A, b)
            α, β, tol, maxit = 1.0, 1e-18, 1e-8, 200
            Atb, AtA = A'b, A'A
            X0 = (AtA + norm(2Atb, Inf) * I) \ Atb
            X, Λ = copy(X0), Diagonal(zeros(size(A, 2)))
            crit, ee = 1.0, 0
            while crit > tol && ee <= maxit
                Λ = Diagonal(α ./ (β .+ abs2.(X0)))
                X = (AtA + Λ) \ Atb
                crit = norm(X - X0, 1) / norm(X0, 1)
                ee += 1
                X0 = X
            end
            return X, inv(AtA + Λ)
        end
        Random.seed!(9)
        for T in (Float64, ComplexF64), (nm, nr) in ((30, 20), (15, 25))
            σ = 0.05
            u = zeros(T, nr); u[[4, 13]] .= (1, -2)
            # dense operator, R = σ² I (τn = 1)
            H = randn(T, nm, nr)
            y = H * u + σ * randn(T, nm)
            s = solve(InverseProblem(H, y), Bayesian(Lq(2.0)), ElementWise(); R = σ^2 * I, covariance = true)
            Xr, Pr = rvr_ref(H / σ, y / σ)
            @test s.u ≈ Xr rtol = 1e-6
            @test s.covariance ≈ Pr rtol = 1e-6
            @test s.hyper.τn == 1
            # block-diagonal operator: one RVR per frequency
            Hs = [randn(T, nm, nr) for _ in 1:3]
            Y = reduce(hcat, [Hj * u + σ * randn(T, nm) for Hj in Hs])
            s = solve(InverseProblem(Hs, Y), Bayesian(Lq(2.0)), ElementWise(); R = σ^2 * I, covariance = true)
            for j in 1:3
                Xr, Pr = rvr_ref(Hs[j] / σ, Y[:, j] / σ)
                @test s.u[:, j] ≈ Xr rtol = 1e-5
                @test s.covariance[j] ≈ Pr rtol = 1e-5
            end
        end
        @test_throws ArgumentError solve(InverseProblem(randn(5, 4), randn(5)), Bayesian(Lq(1.0)), ElementWise())
        @test_throws ArgumentError solve(InverseProblem(randn(5, 4), randn(5)), Bayesian(Lq(2.0); τs = 1.0), ElementWise())
    end

    @testset "Bayesian hyperparameters and credible intervals" begin
        @test_throws ArgumentError Bayesian(Lq(1.0); τn = -1.0)
        @test_throws ArgumentError Bayesian(Lq(1.0); τn = [1.0])
        # Gaussian prior with fixed exact hyperparameters: exact Gaussian posterior
        #   real:    u ~ N(0, 1)  ⇔ exp(-τs u²/2) with τs = 1 ; noise N(0, σ²), τn = 1
        #   complex: u ~ CN(0, 1) ⇔ exp(-τs |u|²) with τs = 1  ; noise CN(0, σ²), τn = 1
        Random.seed!(7)
        n, m, σ = 10, 30, 0.5
        for (T, τs) in ((Float64, 1.0), (ComplexF64, 1.0))
            H = randn(T, m, n)
            reg = Bayesian(Lq(2.0); τn = 1.0, τs = τs, q = Fixed())
            Σexact = inv(H'H / σ^2 + I)
            covered, total = 0, 0
            for _ in 1:200
                u = randn(T, n)
                y = H * u + σ * randn(T, m)
                s = solve(InverseProblem(H, y), reg; R = σ^2 * I, covariance = true)
                @test s.u ≈ Σexact * (H'y) / σ^2
                @test s.covariance ≈ Σexact
                lo, up = credible_interval(s; level = 0.95, nsamples = 4000)
                if T <: Real
                    covered += count(lo .<= u .<= up); total += n
                else
                    covered += count(real(lo) .<= real(u) .<= real(up)) +
                               count(imag(lo) .<= imag(u) .<= imag(up))
                    total += 2n
                end
            end
            @test 0.93 < covered / total < 0.97
            # samples u = μ + S z, Σ = S Sᴴ: empirical mean and covariance
            y = H * randn(T, n) + σ * randn(T, m)
            s = solve(InverseProblem(H, y), reg; R = σ^2 * I, covariance = true)
            X = posterior_samples(s, 100_000; rng = Xoshiro(1))
            @test size(X) == (n, 100_000)
            @test isapprox(vec(sum(X; dims = 2)) / 100_000, s.u; atol = 0.02)
            Xc = X .- s.u
            @test isapprox(Xc * Xc' / 100_000, Σexact; atol = 0.02)
        end
        # estimated vs fixed hyperparameters
        Random.seed!(1)
        H = beam_frf((ωn5 + ωn6) / 2)
        u = zeros(ComplexF64, size(H, 2)); u[17] = 1
        y0 = H * u
        y = y0 .+ 0.05norm(y0) / sqrt(length(y0)) .* randn(ComplexF64, length(y0))
        prob = InverseProblem(H, y)
        se = solve(prob, Bayesian(Lq(1.0)))
        sf = solve(prob, Bayesian(Lq(1.0); τn = se.hyper.τn, q = Fixed()))
        @test sf.hyper.τn == se.hyper.τn
        # default τn: estimated without R, 1 (exact covariance) when R is given
        σ2 = (0.05norm(y0))^2 / length(y0)
        @test se.hyper.τn != 1
        @test solve(prob, Bayesian(Lq(1.0)); R = σ2 * I).hyper.τn == 1
        @test solve(prob, Bayesian(Lq(1.0); τn = Estimated()); R = σ2 * I).hyper.τn != 1
        @test sf.hyper.q == [1.0]
        @test sf.λ ≈ sf.hyper.τs ./ sf.hyper.τn          # complex: ‖r‖² + λ‖u‖_q^q
        @test argmax(abs.(sf.u)) == 17
    end
end
