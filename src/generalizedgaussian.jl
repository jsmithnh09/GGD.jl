"""
    GeneralizedGaussian( μ, α, β )

The "Generalized Normal" or "Generalized Gaussian" with shape parameter 'β',
scale 'α', and location 'μ' has probability density function

```math
f(x; \\mu \\alpha \\beta ) = \\frac{\\beta}{2\\alpha\\Gamma(1/\\beta)}
\\text{exp}\\left( -|x - \\mu|/\\alpha \\right),

\\alpha = \\sigma \\sqrt{ \\frac{\\Gamma(1/\\beta)}{\\Gamma(3/\\beta)} }

```
where `β = 1` incorporates the Laplacian distribution, `β = 2` is the Normal distribution,
and as `β → ∞`, the distribution approaches Uniform on `[μ-α, μ+α]`.

```julia
GeneralizedGaussian(m, a, b)    # Generalized Gaussian with mu, alpha, and beta.
params(d)                       # Get the parameters (m,a,b)
shape(d)                        # returns β shape
scale(d)                        # returns the scale parameter, α
```

External Links
* [Generalized Normal Distribution on Wikipedia](https://en.wikipedia.org/wiki/Generalized_normal_distribution)
* [Reference implementation paper](https://www.researchgate.net/publication/254282790_Simulation_of_the_p-generalized_Gaussian_distribution)

References
[1] Nadarajah, Saralees (September 2005). "A generalized normal distribution".
    Journal of Applied Statistics. 32 (7): 685&ndash, 694.
[2]  Gonzalez-Farias, G., Molina, J. A. D., & Rodríguez-Dagnino, R. M. (2009).
    Efficiency of the approximated shape parameter estimator in the generalized
    Gaussian distribution. IEEE Transactions on Vehicular Technology, 58(8),
    4214-4223.
"""
struct GeneralizedGaussian{T<:Real}
    μ::T
    α::T
    β::T
    function GeneralizedGaussian{T}(μ::T, α::T, β::T) where {T}
        α > zero(α) || error("standard deviation must be positive.")
        β > zero(β) || error("β must be greater than zero.")
        new{T}(μ, α, β)
    end
end

GeneralizedGaussian(μ::T, α::T, β::T) where {T<:Real} = GeneralizedGaussian{T}(μ, α, β)
GeneralizedGaussian(μ::Real, α::Real, β::Real) = GeneralizedGaussian(promote(μ, α, β)...)
GeneralizedGaussian(μ::Integer, σ::Integer, β::Integer) = GeneralizedGaussian(Float64(μ), Float64(σ), Float64(β))

"""
    GeneralizedGaussian(β)

Builds a default distribution with shape `β` and `μ=0, α=1`.
"""
GeneralizedGaussian(β::Real) = GeneralizedGaussian(zero(β), one(β), β)

"""
    GeneralizedGaussian()

Builds the Normal distribution case, where `μ=0, α=√2, β=2`, or N(0, 1).
"""
GeneralizedGaussian() = GeneralizedGaussian(0.0, √2, 2)

### Conversion
function convert(::Type{GeneralizedGaussian{T}}, μ::Real, σ::Real, β::Real) where {T<:Real}
    GeneralizedGaussian(T(μ), T(α), T(β))
end
function convert(::Type{GeneralizedGaussian{T}}, d::GeneralizedGaussian{S}) where {T<:Real,S<:Real}
    GeneralizedGaussian(T(d.μ), T(d.α), T(d.β))
end

### Parameters
params(d::GeneralizedGaussian) = (d.μ, d.α, d.β)
location(d::GeneralizedGaussian) = d.μ
scale(d::GeneralizedGaussian) = d.α
shape(d::GeneralizedGaussian) = d.β

var(d::GeneralizedGaussian) = (d.α^2) * (gamma(3.0 * inv(d.β)) / gamma(inv(d.β)))
std(d::GeneralizedGaussian) = (d.α) * sqrt(gamma(3.0 * inv(d.β)) / gamma(inv(d.β)))

mean(d::GeneralizedGaussian) = d.μ
median(d::GeneralizedGaussian) = d.μ
mode(d::GeneralizedGaussian) = d.μ

### Statistics
skewness(d::GeneralizedGaussian{T}) where {T} = zero(T)
kurtosis(d::GeneralizedGaussian) = gamma(5.0 * inv(d.β)) * gamma(inv(d.β)) / (gamma(3.0 * inv(d.β))^2) - 3.0
entropy(d::GeneralizedGaussian) = inv(d.β) - log(d.β / (2.0 * d.α * gamma(inv(d.β))))

### Evaluation
function pdf(d::GeneralizedGaussian{T}, x::Real) where {T<:Real}
    if x == -Inf || x == Inf
        return zero(T)
    else
        A = d.β / (2.0 * d.α * gamma(1.0 / d.β))
        return A * exp(-(abs(x - d.μ) / d.α)^d.β)
    end
end
logpdf(d::GeneralizedGaussian, x::Real) = log(pdf(d, x))


"""
    cdf(d, x)
Calculates the CDF of the distribution. To determine the CDF, the incomplete
gamma function is required. The CDF  of the Gamma distribution provides this,
with the necessary 1/Γ(a) normalization.
"""
function cdf(d::GeneralizedGaussian{T}, x::Real) where {T<:Real}
    if x == -Inf || x == Inf
        return zero(T)
    else
        v = cdf(Gamma(inv(d.β), 1), (abs(x - μ) / d.α)^d.β) * inv(2)
        return typeof(v)(1 / 2) + sign(x - μ) * v
    end
end


"""
    rand(rng, d)
Extract a sample from the Generalized Gaussian distribution `d`. The sampling
procedure is implemented from from [2].
"""
function rand(rng::AbstractRNG, d::GeneralizedGaussian)
    # utilizing the sampler from the Gamma distribution.
    g = Gamma(inv(d.β), 1)
    # random variable with value -1 or 1 with probability (1/2).
    b = 2.0 * rand(Bernoulli()) - 1
    return d.μ + inv(sqrt(d.α)) * rand(rng, g)^inv(d.β) * b
end

# multi-sample case
function rand(rng::AbstractRNG, d::GeneralizedGaussian, dims::Dims)
    out = Array{eltype(params(d))}(undef, dims)
    @inbounds for i in eachindex(out)
        out[i] = rand(rng, d)
    end
    out
end

# various function signatures for sampling the distribution.
rand(d::GeneralizedGaussian) = rand(GLOBAL_RNG, d)
rand(d::GeneralizedGaussian, dims::Dims) = rand(GLOBAL_RNG, d, dims)
rand(d::GeneralizedGaussian, dims::Int) = rand(GLOBAL_RNG, d, (dims))
rand(rng::AbstractRNG, d::GeneralizedGaussian, dim::Int, dims::Int...) = rand(rng, d, (dim, dims...))
rand(d::GeneralizedGaussian, dims::Int...) = rand(GLOBAL_RNG, d, dims...)


function Zn(X::AbstractArray, β::Real)
    n = size(X, 1)
    S1 = 1 / n * sum(abs.(X) .^ (2β))
    S2 = 1 / n * sum(abs.(X) .^ β)
    return (S1 / (S2^2)) - (β + 1)
end

function Zp(X::AbstractArray, β::Real)
    n = size(X, 1)
    S1 = sum(abs.(X) .^ β) # sum
    S2 = sum(abs.(X) .^ (2β))
    L1 = sum(abs.(X) .^ β .* log.(abs.(X)))    # log
    L2 = sum(abs.(X) .^ (2β) .* log.(abs.(X)))
    num1 = (2 / n * L2) * (1 / n * S1)^2
    denom = (1 / n * S1)^4
    num2 = (1 / n * L1) * (1 / n * S2) * (2 / n * S1)
    return (num1 / denom) - (num2 / denom) - 1
end

"""
    β = gcmsearch(X::AbstractArray, βi::Real)

Globally convergent method for β-parameter estimation using Newton-Raphson
iterative search from [1].

[1] Song, Kai-Sheng. "A globally convergent and consistent method 
    for estimating the shape parameter of a generalized Gaussian distribution." 
    IEEE Transactions on Information Theory 52.2 (2006): 510-527.
"""


function gcmsearch(X::AbstractArray{T}, βi::Real, n::Int) where {T<:AbstractFloat}
    iter = 1
    βcur = βi - Zn(X, βi) / Zp(X, βi)
    while (iter < n) && !isapprox(βcur, βi, atol=√eps(T))
        βi = copy(βcur)
        βcur = βi - (Zn(X, βi) / Zp(X, βi))
        iter += 1
    end
    βcur
end

### using 1e3 trials for "n" trials.
gcmsearch(X, βi) = gcmsearch(X, βi, 1_000)

### using |μ| / σ as an initial estimate when not provided.
gcmsearch(X) = gcmsearch(X, (sum(abs, X) / length(X)) / std(X) + 3, 1_000)

# confidence interval for the estimate of the shape parameter via `gcmsearch`.
# defaulting to z-score for 95 %.
function gcmci(βest::Real, N::Int; z::AbstractFloat=1.96)
    A = coth(acoth(sqrt(βest + 1)) + ((1/(sqrt(2*N)))*z))^2
    B = coth(acoth(sqrt(βest + 1)) - ((1/(sqrt(2*N)))*z))^2
    (A-1, B-1)
end


    
