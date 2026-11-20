using SourceEstimation
using LinearAlgebra
using Random
using SparseArrays
using Test

# Simply supported beam, modal synthesis of the transfer functions
function beam_frf(ω; nm = 30, nu = 25, nmodes = 30, ξ = 0.01)
    xm = range(0, 1, length = nm + 2)[2:end-1]
    xr = range(0, 1, length = nu + 2)[2:end-1]
    ωn = [(n * π)^2 * 50 for n in 1:nmodes]
    return [sum(sin(n * π * xm[i]) * sin(n * π * xr[j]) / (ωn[n]^2 - ω^2 + 2im * ξ * ωn[n] * ω)
                for n in 1:nmodes) for i in 1:nm, j in 1:nu]
end

# Reference inverse by a factorization-based solve (no explicit inverse): A \ I
solve_identity(A) = Matrix(A) \ Matrix{eltype(A)}(I, size(A)...)

const ωn5, ωn6 = (5π)^2 * 50, (6π)^2 * 50

@testset "SourceEstimation.jl" begin
    @testset "Problem construction" begin
        H = randn(10, 6)
        @test size(InverseProblem(H, randn(10)).y) == (10,)
        @test_throws DimensionMismatch InverseProblem(H, randn(9))
        # space–time: H acts on vec(U), U is 3 × 2
        p = InverseProblem(randn(10, 6), randn(5, 2))
        @test SourceEstimation.unknown_size(p) == (3, 2)
        # space–frequency, vector of blocks or 3D array
        Hs = [randn(4, 3) for _ in 1:5]
        @test SourceEstimation.unknown_size(InverseProblem(Hs, randn(4, 5))) == (3, 5)
        @test SourceEstimation.unknown_size(InverseProblem(randn(4, 3, 5), randn(4, 5))) == (3, 5)
    end

    @testset "Regularization terms" begin
        @test_throws ArgumentError Lq(0.0)
        @test_throws ArgumentError Lq(3.0)
        @test_throws ArgumentError Lpq(2, 2.5)
        U = [1.0 -2.0; 0.0 3.0]
        @test SourceEstimation.penalty(Lq(1.0), U) ≈ 6
        @test SourceEstimation.penalty(Lpq(2, 1), U) ≈ sqrt(5) + 3
        @test SourceEstimation.penalty(Lpq(2, 2), U) ≈ SourceEstimation.penalty(Lq(2.0), U)
        G = Groups(1:1 => Lq(1.0), 2:2 => Lq(2.0))
        R = SourceEstimation.Regularization(G, 2)
        @test R.L isa SourceEstimation.IdentitySmoothing
        @test SourceEstimation.penalty(R, U) ≈ 3 + 9
        # one smoothing matrix per group, assembled into a global operator
        D = [1.0 -1.0]
        R = SourceEstimation.Regularization(Groups([1, 2] => Lq(2.0; L = D), [3] => Lq(1.0)), 3)
        @test Matrix(R.L.L) == [1 -1 0; 0 0 1]
        @test R.rows == [[1], [2]]
        @test_throws DimensionMismatch SourceEstimation.groups(Groups([1] => Lq(2.0; L = D), [2] => Lq(1.0)), 2)
        @test_throws ArgumentError SourceEstimation.groups(Groups(1:2 => Lq(1.0), 2:3 => Lq(2.0)), 3)
        @test_throws ArgumentError SourceEstimation.groups(Groups(1:2 => Lq(1.0)), 3)
    end

    @testset "Finite-difference operators" begin
        @test fd_weights(1.0, 3, 1) ≈ [-1.5 2 -0.5; -0.5 0 0.5; 0.5 -2 1.5]
        @test fd_weights(0.5, 3, 2) ≈ fill(4.0, 3) * [1 -2 1]
        @test_throws ArgumentError fd_weights(1.0, 3, 3)
        # 1D: one-sided and centered schemes, exact for polynomials of degree < npts
        @test Matrix(diff_operator(4, 1; npts = 2)) == [-1 1 0 0; 0 -1 1 0; 0 0 -1 1]
        @test diff_operator(5, 0) == sparse(1.0I, 5, 5)
        @test_throws ArgumentError diff_operator(5, 1; npts = 4)
        x = range(0, 2, length = 13)
        h = step(x)
        @test diff_operator(13, 1; h) * x .^ 2 ≈ 2x
        @test diff_operator(13, 2; h) * x .^ 2 ≈ fill(2, 11)            # 3-point: interior points
        @test diff_operator(13, 1; npts = 5, h) * x .^ 4 ≈ 4x .^ 3
        @test size(diff_operator(13, 2; npts = 3)) == (11, 13)
        @test diff_operator(13, 2; npts = 5, h) * x .^ 3 ≈ 6x                 # centered, n × n
        # 2D structured mesh, points in random order
        Random.seed!(10)
        P = shuffle([(xi, yi) for xi in range(0, 0.6, length = 7) for yi in range(0, 0.4, length = 5)])
        C = [first.(P) last.(P)]
        f = C[:, 1] .^ 2 .+ 3C[:, 2] .^ 2
        @test diff_operator(C, 1) * f ≈ 2C[:, 1]
        @test diff_operator(C, 2) * f ≈ 6C[:, 2]
        @test diff_operator(C, 2; order = 2) * f ≈ fill(6, 7 * (5 - 2))
        @test size(diff_operator(C, 1; npts = 2)) == (35 - 5, 35)
        @test_throws ArgumentError diff_operator([0.0 0; 1 0; 3 0], 1)    # non-uniform spacing
        # space–frequency / space–time operators of vectorized problems
        nu, ns = 6, 5
        Ux = randn(ComplexF64, nu, ns)
        Ls, Lx = diff_operator(ns, 1; npts = 2), diff_operator(nu, 2)
        @test space_sample_operator(nu, ns; Lc = Ls) * vec(Ux) ≈ vec(Ux * transpose(Ls))
        @test space_sample_operator(nu, ns; L = Lx) * vec(Ux) ≈ vec(Lx * Ux)
        @test space_sample_operator(nu, ns; L = Lx, Lc = Ls) * vec(Ux) ≈ vec(Lx * Ux * transpose(Ls))
        @test size(space_sample_operator(nu, ns; Lc = Ls)) == ((ns - 1) * nu, ns * nu)
        @test space_sample_operator(nu, ns) == sparse(1.0I, nu * ns, nu * ns)
        @test_throws DimensionMismatch space_sample_operator(nu, ns; Lc = diff_operator(ns + 1, 1))
        @test sample_rows([2, 4], nu, ns) == [2, 4, 8, 10, 14, 16, 20, 22, 26, 28]
        rows = [2, 3, 4]
        @test space_sample_operator(3, ns; Lc = Ls) * vec(Ux)[sample_rows(rows, nu, ns)] ≈
              vec(Ux[rows, :] * transpose(Ls))
        # space–frequency operator of MSSP 2018: Lf of Eq. (11), L = (I_N ⊗ Ls)(Lf ⊗ I_M) of
        # Eq. (15), regularization term ‖L F‖_{p,q}^q of Eq. (14)
        M, N, Δf = 6, 7, 2.0
        Lf11 = zeros(N, N); Lf11[1, 1:3] = [-1.5, 2, -0.5]; Lf11[N, N-2:N] = [0.5, -2, 1.5]
        for j in 2:N-1; Lf11[j, j-1], Lf11[j, j+1] = -0.5, 0.5; end
        Lsp, Lfr = diff_operator(M, 2), diff_operator(N, 1; h = Δf)
        @test Matrix(Lfr) ≈ Lf11 ./ Δf
        L15 = kron(sparse(I, N, N), Lsp) * kron(Lfr, sparse(I, M, M))
        @test L15 ≈ space_sample_operator(M, N; L = Lsp, Lc = Lfr)
        Fm = randn(ComplexF64, M, N)
        Ym = Lsp * Fm * transpose(Lfr)
        @test sum(norm(Ym[i, :], 1.5)^0.7 for i in axes(Ym, 1)) ≈
              SourceEstimation.penalty(Lpq(1.5, 0.7; samples = N), reshape(L15 * vec(Fm), :, 1))
        # use as smoothing operator of a group
        H = randn(40, 35)
        y = H * f
        L = diff_operator(C, 1)
        s = solve(InverseProblem(H, y), Additive(Groups(1:35 => Lq(2.0; L)); λ = 1e-3))
        @test s.u ≈ (H'H + 1e-3 * L'L) \ (H'y)
    end

    @testset "Formulations" begin
        @test_throws ArgumentError Additive(Lq(2.0); λ = -1.0)
        @test_throws ArgumentError Additive(Lq(2.0); λ = Heuristic())
        @test SourceEstimation.default_algorithm(Multiplicative(Lq(1.0))) isa IR
        @test SourceEstimation.default_algorithm(Bayesian(Lq(1.0))) isa MAP
    end

    @testset "Weighted solves" begin
        Random.seed!(0)
        IP = SourceEstimation
        # every regime × smoothing operator × field, against the explicit formulas
        for T in (Float64, ComplexF64), (nm, nu) in ((12, 8), (8, 8), (5, 8)), ns in (1, 3)
            Lsq = diagm(nu, nu, 0 => fill(2.0, nu), 1 => -ones(nu - 1))       # square invertible
            Lrect = diagm(nu, nu, 0 => -ones(nu), 1 => ones(nu - 1))[1:nu-1, :]   # rectangular
            for (Lop, Lm) in ((IP.IdentitySmoothing(), Matrix(1.0I, nu, nu)),
                              (IP.MatrixSmoothing(sparse(Lsq)), Lsq),
                              (IP.MatrixSmoothing(sparse(Lrect)), Lrect))
                Lg = kron(I(ns), Lm)
                # dense operator acting on vec(U) (space–time layout when ns > 1)
                H, y = randn(T, nm * ns, nu * ns), randn(T, nm, ns)
                c = IP.linear_cache(InverseProblem(H, ns == 1 ? vec(y) : y))
                @test (c isa IP.Underdetermined) == (nm < nu)
                D = rand(size(Lm, 1), ns) .+ 0.1
                A = H'H + Lg' * Diagonal(vec(D)) * Lg
                @test vec(IP.weighted_solve(c, Lop, D)) ≈ A \ (H' * vec(y))
                @test IP.inverse_apply(IP.posterior_factor(c, Lop, D), Matrix{T}(I, nu * ns, nu * ns)) ≈ solve_identity(A)
                @test IP.precision_matrix(c, Lop, D, 2.0) ≈ 2A
                # block-diagonal operator
                Hs = [randn(T, nm, nu) for _ in 1:ns]
                cb = IP.linear_cache(InverseProblem(Hs, y))
                U = IP.weighted_solve(cb, Lop, D)
                Fb = IP.posterior_factor(cb, Lop, D)
                Σ = [IP.inverse_apply(Fb[j], Matrix{T}(I, nu, nu)) for j in 1:ns]
                for j in 1:ns
                    Aj = Hs[j]'Hs[j] + Lm' * Diagonal(D[:, j]) * Lm
                    @test U[:, j] ≈ Aj \ (Hs[j]' * y[:, j])
                    @test Σ[j] ≈ solve_identity(Aj)
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
        c = SourceEstimation.linear_cache(InverseProblem(H, y))
        W = rand(12, 1) .+ 0.5
        Lid = SourceEstimation.IdentitySmoothing()
        Lmat = SourceEstimation.MatrixSmoothing(sparse(1.0I, 12, 12))
        # the SVD and generalized SVD paths describe the same problem
        for rule in (LCurve(), GCV(), BayesEstimator())
            λ1 = SourceEstimation.select_parameter(rule, c, Lid, W)
            λ2 = SourceEstimation.select_parameter(rule, c, Lmat, W)
            @test λ1 > 0
            @test isapprox(λ1, λ2; rtol = 1e-6)
        end
        # 1D refinement: Brent's method, and refined grid search better than the grid alone
        x, fx = SourceEstimation.brent(x -> (x - 0.3)^2 + 1, -2, 3)
        @test x ≈ 0.3 atol = 1e-7
        @test fx ≈ 1
        @test SourceEstimation.grid_minimize(x -> (x - 0.3)^2 + 1, range(-2, 3; length = 11)) ≈ 0.3 atol = 1e-7
        sp = SourceEstimation.spectrum(c, Lid, W)
        t = sort(log.(SourceEstimation.lambda_grid(sp.γ)))
        for rule in (LCurve(), GCV(), BayesEstimator())
            J = SourceEstimation.criterion(rule, sp, length(y))
            λ = SourceEstimation.select_parameter(rule, sp, length(y))
            @test J(log(λ) / 2) <= minimum(J.(t))
        end
        s = solve(InverseProblem(H, y), Additive(Lq(2.0); λ = BayesEstimator()))
        @test s.λ > 0 && s.converged
        # non-square L: the null space of L is left unregularized
        Dx = diagm(12, 12, 0 => -ones(12), 1 => ones(11))[1:11, :]
        sp = SourceEstimation.spectrum(c, SourceEstimation.MatrixSmoothing(sparse(Dx)), ones(11, 1))
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
        # posterior covariance: one nu × nu matrix per frequency
        s = solve(prob, Additive(Lq(2.0); λ = 1e-12); precision = true)
        @test length(s.precision) == length(ωs) && length(posterior_covariance(s)) == length(ωs)
        @test size(posterior_std(s)) == size(U)
        @test size(posterior_samples(s, 10)) == (size(U)..., 10)
        lo, up = credible_interval(s; nsamples = 500)
        @test size(lo) == size(U) && all(real(lo) .<= real(s.u) .<= real(up))
    end

    @testset "Operator along the samples Lf (MSSP 2018)" begin
        IP = SourceEstimation
        Random.seed!(14)
        nm, nu, nf = 14, 10, 8
        Hs = [randn(ComplexF64, nm, nu) for _ in 1:nf]
        U = zeros(ComplexF64, nu, nf); U[4, :] .= 1
        Y = reduce(hcat, [Hs[j] * U[:, j] for j in 1:nf]) .+ 0.01randn(ComplexF64, nm, nf)
        Lf = diff_operator(nf, 1)
        pb, pv = InverseProblem(Hs, Y), InverseProblem(blockdiag(sparse.(Hs)...), vec(Y))
        @test IP.sample_coupling(Lpq(2, 0.5; Lc = Lf)) isa IP.SmoothedSamples
        @test IP.sample_coupling(Groups(1:5 => Lq(1.0), 6:10 => Lpq(2, 0.5; Lc = Lf))) isa IP.SmoothedSamples
        @test IP.strategy(pb, Multiplicative(Lpq(2, 0.5; Lc = Lf))) isa IP.Assembled
        # native Lf = vectorized problem with L = kron(Lf, I)
        a = solve(pb, Multiplicative(Lpq(2, 0.5; Lc = Lf)); init = 1e-2, precision = true)
        b = solve(pv, Multiplicative(Lpq(2, 0.5; L = space_sample_operator(nu, nf; Lc = Lf), samples = nf));
                  init = 1e-2, precision = true)
        @test vec(a.u) ≈ b.u && a.λ ≈ b.λ && size(a.u) == (nu, nf)
        # sparse precision; marginal blocks (block solves) or full covariance
        @test a.precision isa SparseMatrixCSC && a.precision ≈ b.precision
        Σb = posterior_covariance(b)                                     # vectorized: one block
        @test posterior_covariance(a; full = true) ≈ Σb ≈ solve_identity(Matrix(a.precision))
        Σa = posterior_covariance(a)
        @test length(Σa) == nf && size(Σa[1]) == (nu, nu)
        @test all(Σa[j] ≈ Σb[(j-1)*nu+1:j*nu, (j-1)*nu+1:j*nu] for j in 1:nf)
        @test posterior_std(a) ≈ reshape(sqrt.(real(diag(Σb))), nu, nf)
        lo, up = credible_interval(a; nsamples = 200)
        @test size(lo) == (nu, nf)
        # joint samples (sparse Cholesky R⁻¹ z): empirical covariance → full covariance
        X = reshape(posterior_samples(a, 40_000; rng = Xoshiro(2)), nu * nf, :)
        Xc = X .- vec(a.u)
        @test isapprox(Xc * Xc' / size(X, 2), Σb; rtol = 0.05)
        # space–time (dense): marginal blocks per time step from the dense covariance
        Hst, Yst = randn(nm * 4, nu * 4), randn(nm, 4)
        sm = solve(InverseProblem(Hst, Yst), Additive(Lpq(2, 1); λ = 0.1); precision = true)
        Σm, Σf = posterior_covariance(sm), posterior_covariance(sm; full = true)
        @test length(Σm) == 4 && Σm[3] ≈ Σf[2nu+1:3nu, 2nu+1:3nu] && Σf ≈ solve_identity(Matrix(sm.precision))
        # groups, per-point noise variances
        r = 0.01 .+ rand(nm)
        a = solve(pb, Multiplicative(Groups(1:5 => Lpq(2, 0.5; Lc = Lf), 6:10 => Lq(1.0))); R = r, init = 1e-2)
        b = solve(pv, Multiplicative(Groups(sample_rows(1:5, nu, nf) => Lpq(2, 0.5; L = space_sample_operator(5, nf; Lc = Lf), samples = nf),
                                            sample_rows(6:10, nu, nf) => Lq(1.0))); R = repeat(r, nf), init = 1e-2)
        @test vec(a.u) ≈ b.u
        # space–time problem (dense convolution matrix), spatial and temporal operators
        nt = 6
        Ht, Yt = randn(nm * nt, nu * nt), randn(nm, nt)
        Lsp, Lt = diff_operator(nu, 2; npts = 5), diff_operator(nt, 1)
        a = solve(InverseProblem(Ht, Yt), Additive(Lpq(2, 1; L = Lsp, Lc = Lt); λ = 0.1))
        b = solve(InverseProblem(Ht, vec(Yt)),
                  Additive(Lpq(2, 1; L = space_sample_operator(nu, nt; L = Lsp, Lc = Lt), samples = nt); λ = 0.1))
        @test vec(a.u) ≈ b.u
        # errors
        @test_throws ArgumentError solve(InverseProblem(Hs[1], Y[:, 1]), Multiplicative(Lpq(2, 0.5; Lc = Lf)))
        @test_throws ArgumentError Lpq(2, 0.5; Lc = Lf, samples = 3)
        @test_throws DimensionMismatch solve(pb, Multiplicative(Lpq(2, 0.5; Lc = diff_operator(nf + 1, 1))))
        @test_throws ArgumentError solve(pv, Multiplicative(Lpq(2, 0.5; Lc = Lf)))
        # sparse heuristic on H L⁻¹ (invertible L) = dense heuristic (generalized SVD)
        Hov = randn(ComplexF64, 30, 12)
        S = IP.MatrixSmoothing(sparse(diff_operator(12, 1) + 0.5I))
        cs = IP.linear_cache(InverseProblem(sparse(Hov), randn(ComplexF64, 30)))
        cd = IP.linear_cache(InverseProblem(Hov, randn(ComplexF64, 30)))
        @test IP.select_parameter(Heuristic(), cs, S, ones(12, 1)) ≈ IP.select_parameter(Heuristic(), cd, S, ones(12, 1)) rtol = 1e-3
    end

    @testset "Multiplicative ℓp,q regularization (port of lpq_reg)" begin
        IP = SourceEstimation
        # direct port of lpq_reg, with the IR weights and the damping rule of the package
        function lpq_reg_ref(A, b, p, q, r)
            S = Diagonal(1 ./ sqrt.(r))
            A, b = [S * a for a in A], S * b
            nf = length(A); nb, nx = size(A[1])
            Af = cat(A...; dims = (1, 2)); bf = b[:]
            AtA, Atb = Af'Af, Af'bf
            λ = norm(2Atb, Inf)
            X0 = (AtA + λ * I) \ Atb
            ε = IP.damping(X0)
            W = zeros(nx, nf)
            X = similar(X0)
            for _ in 1:50
                s = IP.reweight!(W, Lpq(p, q), reshape(X0, nx, nf), ε)
                λ = norm(bf - Af * X0)^2 / s
                X = (AtA + λ * Diagonal(vec(W))) \ Atb
                crit = norm(X - X0, 1) / norm(X0, 1)
                X0 = X
                crit < 1e-8 && break
            end
            P = [solve_identity(A[f]'A[f] + λ * Diagonal(W[:, f])) for f in 1:nf]
            return reshape(X, nx, nf), P, λ
        end
        Random.seed!(13)
        ωs = range(2000, 9000, length = 5)
        Hs = [beam_frf(ω) for ω in ωs]
        U = zeros(ComplexF64, 25, length(ωs)); U[17, :] .= 1
        Y0 = reduce(hcat, [Hs[j] * U[:, j] for j in eachindex(ωs)])
        r = (0.05 * maximum(abs, Y0))^2 .* (0.5 .+ rand(size(Y0, 1)))   # one variance per point
        Y = Y0 .+ sqrt.(r) .* randn(ComplexF64, size(Y0))
        Xr, Pr, λr = lpq_reg_ref(Hs, Y, 2.0, 0.5, r)
        s = solve(InverseProblem(Hs, Y), Multiplicative(Lpq(2, 0.5)); R = r, init = LambdaMax(),
                  maxiter = 50, tol = 1e-8, precision = true)   # stopping rule of the reference (crit < 1e-8)
        @test s.u ≈ Xr rtol = 1e-8
        @test s.λ ≈ λr rtol = 1e-8
        @test all(posterior_covariance(s)[f] ≈ Pr[f] for f in eachindex(ωs))
        # per-point variances: vector, Diagonal or diagonal matrix nm × nm
        for R in (Diagonal(r), Matrix(Diagonal(r)))
            @test solve(InverseProblem(Hs, Y), Multiplicative(Lpq(2, 0.5)); R, init = LambdaMax(), maxiter = 50, tol = 1e-8).u ≈ s.u
        end
        @test_throws ArgumentError Additive(Lq(1.0); λ = LambdaMax())
    end

    @testset "Posterior precision and sampling (Cholesky / perturbation–optimization)" begin
        IP = SourceEstimation
        Random.seed!(15)
        # over-determined problems (dense Cholesky): standard sampling u + R⁻¹ z / √τn
        for T in (Float64, ComplexF64)
            H, y = randn(T, 12, 6), randn(T, 12)
            s = solve(InverseProblem(H, y), Additive(Lq(2.0); λ = 0.5); precision = true)
            @test s.posterior.F isa LinearAlgebra.Cholesky
            Σ = solve_identity(Matrix(s.precision))
            X = posterior_samples(s, 60_000; rng = Xoshiro(4))
            Xc = X .- s.u
            @test isapprox(Xc * Xc' / size(X, 2), Σ; rtol = 0.05)
        end
        # under-determined problems (Woodbury identity), with L = I and an invertible L
        for T in (Float64, ComplexF64), Lt in (Lq(2.0), Lq(2.0; L = diff_operator(12, 1) + 0.5I))
            H, y = randn(T, 6, 12), randn(T, 6)
            s = solve(InverseProblem(H, y), Additive(Lt; λ = 0.5); precision = true)
            @test s.posterior.F isa IP.WoodburyFactor
            Σ = solve_identity(Matrix(s.precision))
            @test posterior_covariance(s) ≈ Σ
            @test posterior_std(s) ≈ sqrt.(real(diag(Σ)))
            X = posterior_samples(s, 60_000; rng = Xoshiro(3))
            Xc = X .- s.u
            @test isapprox(Xc * Xc' / size(X, 2), Σ; rtol = 0.05)
        end
        @test_throws ArgumentError posterior_std(solve(InverseProblem(randn(6, 4), randn(6)), Additive(Lq(2.0); λ = 0.5)))
    end

    @testset "Quadratic additive regularization: closed form" begin
        Random.seed!(16)
        H, y = randn(20, 12), randn(20)
        prob = InverseProblem(H, y)
        # fixed λ, multi-parameter λ and parameter choice: one linear solve, no IR iteration
        s = solve(prob, Additive(Lq(2.0); λ = 0.3))
        @test s.iterations == 0 && s.converged && s.u ≈ (H'H + 0.3I) \ (H'y)
        s = solve(prob, Additive(Groups(1:6 => Lq(2.0), 7:12 => Lq(2.0)); λ = [0.1, 1.0]))
        @test s.iterations == 0 && s.u ≈ (H'H + Diagonal([fill(0.1, 6); fill(1.0, 6)])) \ (H'y)
        s = solve(prob, Additive(Lq(2.0); λ = GCV()))
        @test s.iterations == 0 && s.u ≈ (H'H + s.λ * I) \ (H'y)
        # not quadratic: IR iterations
        @test solve(prob, Additive(Lq(1.0); λ = 0.3)).iterations > 0
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
        s = solve(prob, Additive(Lq(2.0); λ = 0.3); R = v, precision = true)
        @test posterior_covariance(s) ≈ solve_identity(H' * (H ./ v) + 0.3I)
        @test s.precision ≈ H' * (H ./ v) + 0.3I
        @test posterior_std(s) ≈ sqrt.(diag(solve_identity(H' * (H ./ v) + 0.3I)))
        @test_throws ArgumentError posterior_std(solve(prob, Additive(Lq(2.0); λ = 0.3)))
    end

    @testset "Generalized Gaussian priors (real / complex)" begin
        IP = SourceEstimation
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
            return X, solve_identity(AtA + Λ)
        end
        Random.seed!(9)
        for T in (Float64, ComplexF64), (nm, nu) in ((30, 20), (15, 25))
            σ = 0.05
            u = zeros(T, nu); u[[4, 13]] .= (1, -2)
            # dense operator, R = σ² I (τn = 1)
            H = randn(T, nm, nu)
            y = H * u + σ * randn(T, nm)
            s = solve(InverseProblem(H, y), Bayesian(Lq(2.0)), ElementWise(); R = σ^2 * I, precision = true)
            Xr, Pr = rvr_ref(H / σ, y / σ)
            @test s.u ≈ Xr rtol = 1e-6
            @test posterior_covariance(s) ≈ Pr rtol = 1e-6
            @test s.hyper.τn == 1
            # block-diagonal operator: one RVR per frequency
            Hs = [randn(T, nm, nu) for _ in 1:3]
            Y = reduce(hcat, [Hj * u + σ * randn(T, nm) for Hj in Hs])
            s = solve(InverseProblem(Hs, Y), Bayesian(Lq(2.0)), ElementWise(); R = σ^2 * I, precision = true)
            for j in 1:3
                Xr, Pr = rvr_ref(Hs[j] / σ, Y[:, j] / σ)
                @test s.u[:, j] ≈ Xr rtol = 1e-5
                @test posterior_covariance(s)[j] ≈ Pr rtol = 1e-5
            end
        end
        @test_throws ArgumentError solve(InverseProblem(randn(5, 4), randn(5)), Bayesian(Lq(1.0)), ElementWise())
        @test_throws ArgumentError solve(InverseProblem(randn(5, 4), randn(5)), Bayesian(Lq(2.0); τs = 1.0), ElementWise())
    end

    @testset "Frequency-by-frequency resolution" begin
        IP = SourceEstimation
        @test IP.sample_coupling(Lq(1.0)) isa IP.IndependentSamples
        @test IP.sample_coupling(Lpq(2, 1)) isa IP.CoupledSamples
        @test IP.sample_coupling(Groups(1:2 => Lq(1.0), 3:4 => Lq(2.0))) isa IP.IndependentSamples
        @test IP.sample_coupling(Groups(1:2 => Lq(1.0), 3:4 => Lpq(2, 1))) isa IP.CoupledSamples
        Random.seed!(11)
        ωs = range(2000, 9000, length = 6)
        Hs = [beam_frf(ω) for ω in ωs]
        U = zeros(ComplexF64, 25, length(ωs)); U[17, :] .= 1
        Y0 = reduce(hcat, [Hs[j] * U[:, j] for j in eachindex(ωs)])
        σ2 = (0.05 .* abs.(Y0)) .^ 2 .+ 1e-30
        Y = Y0 .+ sqrt.(σ2) .* randn(ComplexF64, size(Y0))
        prob = InverseProblem(Hs, Y)
        loop(f, alg; kw...) = [solve(InverseProblem(Hs[j], Y[:, j]), f, alg; R = σ2[:, j], kw...) for j in eachindex(ωs)]
        # Lq: one problem per frequency, one adaptive parameter per frequency
        f = Multiplicative(Lq(0.5))
        s = solve(prob, f; R = σ2)
        ref = loop(f, IR())
        @test s.u ≈ reduce(hcat, [r.u for r in ref])
        @test s.λ ≈ [r.λ for r in ref] && length(s.λ) == length(ωs)
        @test s.hyper.iterations == [r.iterations for r in ref]
        @test s.converged == all(r.converged for r in ref)
        # MPMR with groups of Lq: one parameter per group and per frequency
        f = Multiplicative(Groups(1:12 => Lq(1.0), 13:25 => Lq(0.5)); multiparameter = true)
        @test size(solve(prob, f; R = σ2).λ) == (2, length(ωs))
        # Bayesian MAP: hyperparameters per frequency, covariance per frequency
        f = Bayesian(Lq(1.0))
        s = solve(prob, f; R = σ2, precision = true)
        ref = loop(f, MAP(); precision = true)
        @test s.u ≈ reduce(hcat, [r.u for r in ref])
        @test s.hyper.τn == [r.hyper.τn for r in ref]          # = 1: R exact
        @test size(s.hyper.q) == (1, length(ωs))
        @test length(posterior_covariance(s)) == length(ωs) && posterior_covariance(s)[2] ≈ posterior_covariance(ref[2])
        lo, up = credible_interval(s; nsamples = 200)
        @test size(lo) == size(U)
        # element-wise Bayesian: exactly one RVR per frequency
        s = solve(prob, Bayesian(Lq(2.0)), ElementWise(); R = σ2)
        ref = loop(Bayesian(Lq(2.0)), ElementWise())
        @test s.u ≈ reduce(hcat, [r.u for r in ref])
        @test size(s.λ) == size(U) && size(s.hyper.τs) == size(U)
        # mixed norm: the space–frequency problem is solved as a whole (one parameter)
        s = solve(prob, Multiplicative(Lpq(2, 0.5)); R = σ2)
        @test s.λ isa Real && !haskey(s.hyper, :iterations)
        s = solve(prob, Multiplicative(Groups(1:12 => Lq(1.0), 13:25 => Lpq(2, 0.5))); R = σ2)
        @test s.λ isa Real
    end

    @testset "Sparse operators (vectorized space–frequency problems)" begin
        Random.seed!(12)
        nm, nu, nf = 12, 10, 8
        Hs = [randn(ComplexF64, nm, nu) for _ in 1:nf]
        U = zeros(ComplexF64, nu, nf); U[4, :] .= 1
        Y = reduce(hcat, [Hs[j] * U[:, j] for j in 1:nf]) .+ 0.01randn(ComplexF64, nm, nf)
        Hb = blockdiag(sparse.(Hs)...)
        y = vec(Y)
        v = rand(length(y)) .+ 0.5
        ps, pd = InverseProblem(Hb, y), InverseProblem(Matrix(Hb), y)
        @test ps.H isa SourceEstimation.SparseOperator
        # frequency smoothing L = Lf ⊗ I couples the frequencies
        Lf = kron(diff_operator(nf, 1), sparse(1.0I, nu, nu))
        Hd = Matrix(Hb)
        s = solve(ps, Additive(Lq(2.0; L = Lf); λ = 0.1); R = v, precision = true)
        A = Hd' * (Hd ./ v) + 0.1 * Matrix(Lf' * Lf)
        @test s.u ≈ A \ (Hd' * (y ./ v))
        @test posterior_covariance(s) ≈ solve_identity(A)
        # same results as the dense path
        for (f, alg, kw) in ((Multiplicative(Lq(0.5)), IR(), (init = 1e-2,)),
                             (Bayesian(Lq(1.0)), MAP(), (init = 1e-2,)),
                             (Bayesian(Lq(2.0)), ElementWise(), (;)))
            a, b = solve(ps, f, alg; kw..., precision = true), solve(pd, f, alg; kw..., precision = true)
            @test a.u ≈ b.u rtol = 1e-8
            @test posterior_covariance(a) ≈ posterior_covariance(b) rtol = 1e-6
        end
        # default initializations: heuristic estimate for the OMR and the MAP
        s = solve(ps, Multiplicative(Lq(0.5)))
        @test s.λ > 0 && all(argmax(abs.(reshape(s.u, nu, nf)[:, j])) == 4 for j in 1:nf)
        @test solve(ps, Bayesian(Lq(1.0))).hyper.τn > 0
        # 2-norm condition number estimate, heuristic consistent with the dense one
        A = sparse(Hermitian(Matrix(Hb' * Hb) + 0.1I))
        @test SourceEstimation.condest(A) ≈ cond(Matrix(A)) rtol = 1e-3
        @test SourceEstimation.condest(sparse([1.0 1; 1 1])) == Inf
        Hover = sparse(randn(ComplexF64, 30, 12))
        cs = SourceEstimation.linear_cache(InverseProblem(Hover, randn(ComplexF64, 30)))
        cd = SourceEstimation.linear_cache(InverseProblem(Matrix(Hover), randn(ComplexF64, 30)))
        W = ones(12, 1)
        @test SourceEstimation.select_parameter(Heuristic(), cs, SourceEstimation.IdentitySmoothing(), W) ≈
              SourceEstimation.select_parameter(Heuristic(), cd, SourceEstimation.IdentitySmoothing(), W) rtol = 1e-2
        # mixed norm on the vectorized problem: same model as the space–frequency problem
        pb = InverseProblem(Hs, Y)
        a = solve(pb, Multiplicative(Lpq(2, 0.5)); init = 1e-2)
        b = solve(ps, Multiplicative(Lpq(2, 0.5; samples = nf)); init = 1e-2)
        @test vec(a.u) ≈ b.u && a.λ ≈ b.λ && a.iterations == b.iterations
        @test SourceEstimation.penalty(Lpq(2, 0.5; samples = nf), reshape(vec(U), :, 1)) ≈
              SourceEstimation.penalty(Lpq(2, 0.5), U)
        @test_throws ArgumentError solve(ps, Multiplicative(Lpq(2, 0.5)); init = 1e-2)
        # groups of a vectorized problem, mixed norm combined with a frequency smoothing
        Ls = diff_operator(nf, 1)
        f = Multiplicative(Groups(sample_rows(1:5, nu, nf) => Lpq(2, 0.5; samples = nf),
                                  sample_rows(6:10, nu, nf) => Lpq(2, 1; L = space_sample_operator(5, nf; Lc = Ls),
                                                                   samples = size(Ls, 1))))
        @test solve(ps, f).converged
        @test_throws DimensionMismatch solve(ps, Multiplicative(Groups(sample_rows(1:5, nu, nf) => Lq(1.0),
            sample_rows(6:10, nu, nf) => Lpq(2, 1; L = space_sample_operator(5, nf; Lc = Ls), samples = 3))))
        # SVD-based parameter choice rules are not available
        @test_throws ArgumentError solve(ps, Additive(Lq(1.0); λ = LCurve()))
        @test_throws ArgumentError solve(ps, Multiplicative(Lq(1.0)); init = GCV())
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
            Σexact = solve_identity(H'H / σ^2 + I)
            covered, total = 0, 0
            for _ in 1:200
                u = randn(T, n)
                y = H * u + σ * randn(T, m)
                s = solve(InverseProblem(H, y), reg; R = σ^2 * I, precision = true)
                @test s.u ≈ Σexact * (H'y) / σ^2
                @test posterior_covariance(s) ≈ Σexact
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
            s = solve(InverseProblem(H, y), reg; R = σ^2 * I, precision = true)
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

    @testset "Performance indicators" begin
        uref = [1.0, -2.0, 0.0, 3.0]
        u = [1.0, -1.0, 0.5, 3.0]
        @test relative_error(uref, u) ≈ 100 * 1.5 / 6
        @test relative_error(uref, u; p = 2) ≈ 100 * norm(uref - u) / norm(uref)
        @test correlation(uref, -2uref) ≈ 100
        @test correlation(uref, [0.0, 0.0, 1.0, 0.0]) ≈ 0 atol = 1e-12
        @test_throws DimensionMismatch relative_error(uref, u[1:3])

        # covariances: one per sample (marginal blocks, Kalman filter) or one for all the unknowns
        Uref, U = randn(3, 4), randn(3, 4)
        Σs = [Matrix(Diagonal(fill(Float64(k), 3))) for k in 1:4]
        Σfull = cat(Σs...; dims = (1, 2))
        @test precision_metric(Σs; covariance = true) ≈ sqrt(3 * sum(1:4) / 12)
        @test precision_metric(Σfull; covariance = true) ≈ precision_metric(Σs; covariance = true)
        mse = (sum(abs2, Uref - U) + 3 * sum(1:4)) / 12
        @test accuracy_metric(Uref, U, Σs; covariance = true) ≈ sqrt(mse)
        @test accuracy_metric(Uref, U, Σfull; covariance = true) ≈ sqrt(mse)
        @test accuracy_metric(Uref, U, zeros(12, 12); covariance = true) ≈ sqrt(sum(abs2, Uref - U) / 12)
        @test_throws DimensionMismatch accuracy_metric(Uref, U, Σs[1:3]; covariance = true)

        # precisions (default): dense, Hermitian, sparse, diagonal, blocks — same values as Q⁻¹
        for T in (Float64, ComplexF64)
            B = randn(T, 12, 12)
            Q = B' * B + 12I
            Σ = solve_identity(Q)
            P = precision_metric(Σ; covariance = true)
            @test precision_metric(Q) ≈ P
            @test precision_metric(Hermitian(Q)) ≈ P
            @test precision_metric(sparse(Q)) ≈ P
            @test SourceEstimation.inverse_diagonal(sparse(Q); chunk = 5) ≈ real(diag(Σ))
            @test accuracy_metric(vec(Uref), vec(U), Q) ≈ accuracy_metric(vec(Uref), vec(U), Σ; covariance = true)
            # wrappers storing a single triangle
            @test precision_metric(Hermitian(tril(Q), :L)) ≈ P
            @test precision_metric(Hermitian(sparse(triu(Q)), :U)) ≈ P
            @test precision_metric(Hermitian(sparse(tril(Q)), :L)) ≈ P
            # one block per sample, real or complex
            Qs = [(C = randn(T, 3, 3); C' * C + 3I) for _ in 1:4]
            Σs_T = solve_identity.(Qs)
            @test precision_metric(Qs) ≈ precision_metric(Σs_T; covariance = true)
            @test precision_metric(Qs) ≈ precision_metric(cat(Qs...; dims = (1, 2)))
            Ur, Uc = randn(T, 3, 4), randn(T, 3, 4)
            @test accuracy_metric(Ur, Uc, Qs) ≈ accuracy_metric(Ur, Uc, Σs_T; covariance = true)
            @test precision_metric(Diagonal(T[4, 1])) ≈ sqrt((0.25 + 1) / 2)
        end
        @test precision_metric(Diagonal([4.0, 1.0])) ≈ sqrt((0.25 + 1) / 2)
        @test precision_metric([inv(S) for S in Σs]) ≈ precision_metric(Σs; covariance = true)

        # from a solution: variances by block solves, same values as its precision matrix
        Random.seed!(21)
        for T in (Float64, ComplexF64)
            H, y = randn(T, 12, 6), randn(T, 12)
            s = solve(InverseProblem(H, y), Additive(Lq(2.0); λ = 0.5); precision = true)
            ur = randn(T, 6)
            @test precision_metric(s) ≈ precision_metric(s.precision)
            @test precision_metric(s) ≈ precision_metric(posterior_covariance(s); covariance = true)
            @test accuracy_metric(ur, s) ≈ accuracy_metric(ur, s.u, s.precision)
            @test relative_error(ur, s) == relative_error(ur, s.u)
        end

        ind = [2.0 4.0; 1.0 8.0; 4.0 2.0]
        @test performance_ratio(ind) ≈ [sqrt(2 * 2), sqrt(1 * 4), sqrt(4 * 1)]
        @test_throws ArgumentError performance_ratio([0.0 1.0])
    end

    @testset "Ill-posedness indicators" begin
        rng = Xoshiro(31)
        m, n = 60, 40
        V = Matrix(qr(randn(rng, m, m)).Q)[:, 1:n]
        W = Matrix(qr(randn(rng, n, n)).Q)
        svd_operator(σ) = V * Diagonal(σ) * W'

        # polynomial decay σᵢ = i⁻², smooth field wᵢᴴ u = σᵢ: coefficients σᵢ², Picard satisfied
        σ = (1:n) .^ -2.0
        H = svd_operator(σ)
        u = W * σ
        η = 1e-4
        y = H * u + η * randn(rng, m)
        ip = ill_posedness(H, y; R = η^2 * I)
        @test ip.σ ≈ σ ./ η                                    # whitened operator H / η
        @test !ip.noise_estimated && ip.noise_level == 1
        @test ip.picard_satisfied
        @test 1.2 < ip.picard_slope < 2.5                      # ≈ 2, flattened near the noise floor
        @test 5 <= ip.picard_index <= 10                       # σᵢ² ≈ 2η up to i ≈ 8
        @test ip.decay_rate ≈ 2 atol = 1e-8
        @test ip.severity == :moderate
        @test ip.numerical_rank == n && ip.nullity == 0
        @test ip.condition_number ≈ n^2                        # invariant to the whitening
        @test ill_posedness(H, y).σ ≈ σ                        # R unknown: no whitening
        # unknown noise level: estimated from the residual outside the range of H (m > n)
        ipe = ill_posedness(H, y)
        @test ipe.noise_estimated
        @test 0.7η < ipe.noise_level < 1.4η
        @test abs(ipe.picard_index - ip.picard_index) <= 2
        @test ill_posedness(InverseProblem(H, y); R = η^2 * I).picard_index == ip.picard_index

        # rough field wᵢᴴ u = 1/σᵢ: flat exact coefficients, Picard condition violated
        y = H * (W * (1 ./ σ)) + 1e-10 * randn(rng, m)
        ipv = ill_posedness(H, y; R = 1e-20 * I)
        @test !ipv.picard_satisfied
        @test abs(ipv.picard_slope) < 0.2
        @test ipv.picard_index == n

        # exponential decay: severely ill-posed; complex data
        Hc = svd_operator(exp.(-0.5 .* (1:n))) .* cis(0.3)
        ips = ill_posedness(Hc, Hc * randn(rng, ComplexF64, n))
        @test ips.severity == :severe

        # underdetermined and rank-deficient operators
        ipu = ill_posedness(randn(rng, 10, 20), randn(rng, 10))
        @test length(ipu.σ) == 10 && ipu.numerical_rank == 10 && ipu.nullity == 10
        Hd = randn(rng, 12, 6); Hd[:, 6] .= 0
        ipd = ill_posedness(Hd, randn(rng, 12))
        @test ipd.numerical_rank == 5 && ipd.nullity == 1 && ipd.condition_number == Inf

        # space–frequency problem: one analysis per frequency, per-point noise variances
        Hs = [randn(rng, ComplexF64, 8, 5) for _ in 1:3]
        Y = randn(rng, ComplexF64, 8, 3)
        ipf = ill_posedness(InverseProblem(Hs, Y); R = fill(0.5, 8))
        @test ipf isa Vector{IllPosedness} && length(ipf) == 3
        @test ipf[2].σ ≈ svdvals(Hs[2] ./ sqrt(0.5))

        @test_throws ArgumentError ill_posedness(H, y; window = 4)
        @test sprint(show, MIME("text/plain"), ip) isa String
    end
end

# Makie extension ---------------------------------------------------------------------------------
@testset "Plots without Makie: error hint" begin
    ip = ill_posedness(randn(Xoshiro(1), 20, 10), randn(Xoshiro(2), 20))
    err = try picardplot(ip); nothing catch e; e end
    @test err isa MethodError
    @test occursin("using CairoMakie", sprint(showerror, err))
end

using Makie

@testset "Makie extension: Picard plot" begin
    @test Base.get_extension(SourceEstimation, :SourceEstimationMakieExt) !== nothing
    rng = Xoshiro(3)
    H, y = randn(rng, 30, 12), randn(rng, 30)
    ip = ill_posedness(H, y)
    @test picardplot(ip) isa Makie.Figure
    @test picardplot(H, y; legend = false, figure = (; size = (400, 300))) isa Makie.Figure
    @test picardplot(InverseProblem(H, y); R = 0.5I, window = 3, components = :all, smoothed = false) isa Makie.Figure

    ax = Axis(Figure()[1, 1])
    @test picardplot!(ax, ip; components = 5) === ax
    @test length(ax.scene.plots) == 4 + (ip.picard_index > 0)   # σ, β, ratio, noise level, Picard index
    ax2 = Axis(Figure()[1, 1])
    picardplot!(ax2, ip; noise = false, index = false)
    @test length(ax2.scene.plots) == 3
    @test_throws ArgumentError picardplot(ip; components = :first)

    # zero coefficients (rank-deficient operator, no data) are skipped on the log scale
    @test picardplot(ill_posedness([1.0 0; 0 0; 0 0], zeros(3)); components = :all) isa Makie.Figure

    # space–frequency problem: one analysis per sample
    prob = InverseProblem([randn(rng, ComplexF64, 8, 5) for _ in 1:3], randn(rng, ComplexF64, 8, 3))
    @test picardplot(prob; sample = 2) isa Makie.Figure
    @test_throws ArgumentError picardplot(prob)
end
