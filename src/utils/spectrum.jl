# Spectral data of the standard-form problem ---------------------------------------

"""
    spectrum(cache, L, W)

Spectral description of `min ‖y - Hu‖² + λ ‖W^{1/2} L u‖²` used by the parameter choice rules:

* `γ`: (generalized) singular values of the regularized components;
* `β2`: squared Fourier coefficients `|Uᴴ y|²` of these components;
* `δ0`: part of `‖y‖²` outside the range of `H`;
* `nfree`: number of unregularized components (null space of `L`).

SVD of `H W^{-1/2}` when `L = I`, generalized SVD of `(H, W^{1/2} L)` otherwise.
"""
function spectrum(c::DenseCache, ::IdentitySmoothing, W::AbstractMatrix{<:Real})
    F = svd(c.H ./ transpose(sqrt.(vec(W))))
    β = F.U' * c.y
    return (γ = F.S, β2 = abs2.(β), δ0 = max(sum(abs2, c.y) - sum(abs2, β), 0.0), nfree = 0)
end

function spectrum(c::DenseCache{T}, S::MatrixSmoothing, W::AbstractMatrix{<:Real}) where {T}
    B = Matrix{T}(Diagonal(sqrt.(vec(W))) * global_operator(c, S))
    F = svd(c.H, B)
    ncomp = min(size(c.H, 1), F.k + F.l)
    β = view(F.U, :, 1:ncomp)' * c.y
    reg = (F.k + 1):ncomp
    return (γ = F.alpha[reg] ./ F.beta[reg], β2 = abs2.(β[reg]),
            δ0 = max(sum(abs2, c.y) - sum(abs2, β), 0.0), nfree = F.k)
end

function spectrum(c::BlockCache, S::SmoothingOperator, W::AbstractMatrix{<:Real})
    parts = [spectrum(c.blocks[j], S, view(W, :, j:j)) for j in 1:c.ns]
    return (γ = reduce(vcat, [p.γ for p in parts]), β2 = reduce(vcat, [p.β2 for p in parts]),
            δ0 = sum(p.δ0 for p in parts), nfree = sum(p.nfree for p in parts))
end
