"""
    GeneralizedGaussian( μ, α, β )

The "Generalized Normal" or "Generalized Gaussian" with shape parameter 'β',
scale 'α', and location 'μ' has probability density function

```math
f(x; \\mu \\alpha \\beta ) = \\frac{\\beta}{2\\alpha\\Gamma(1/\\beta)}
\\text{exp}\\left( -|x - \\mu|/\\alpha \\right),

\\alpha = \\sigma \\sqrt{ \\frac{\\Gamma(1/\\beta)}{\\Gamma(3/\\beta)} }

```
where `p = 1` incorporates the Laplacian distribution, `p = 2` is the Normal distribution,
and as `p → ∞`, the distribution approaches Uniform on `[μ-α, μ+α]`.

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
    function GeneralizedGaussian{T}(μ::T, α::T, β::T) where T
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
function convert(::Type{GeneralizedGaussian{T}}, μ::Real, σ::Real, β::Real) where T <: Real
    GeneralizedGaussian(T(μ), T(α), T(β))
end
function convert(::Type{GeneralizedGaussian{T}}, d::GeneralizedGaussian{S}) where {T <: Real, S <: Real}
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
skewness(d::GeneralizedGaussian{T}) where T = zero(T)
kurtosis(d::GeneralizedGaussian) = gamma(5.0 * inv(d.β)) * gamma(inv(d.β)) / (gamma(3.0 * inv(d.β))^2) - 3.0
entropy(d::GeneralizedGaussian) = inv(d.β) - log( d.β / (2.0 * d.α * gamma(inv(d.β))))

### Evaluation
function pdf(d::GeneralizedGaussian{T}, x::Real) where T<:Real
    if x == -Inf || x == Inf
        return zero(T)
    else
        A = d.β / (2.0 * d.α * gamma(1. / d.β) )
        return A * exp( -(abs(x - d.μ) / d.α)^d.β )
    end
end
logpdf(d::GeneralizedGaussian, x::Real) = log(pdf(d, x))


"""
    cdf(d, x)
Calculates the CDF of the distribution. To determine the CDF, the incomplete
gamma function is required. The CDF  of the Gamma distribution provides this,
with the necessary 1/Γ(a) normalization.
"""
function cdf(d::GeneralizedGaussian{T}, x::Real) where T<:Real
    if x == -Inf || x == Inf
        return zero(T)
    else
        v = cdf(Gamma(inv(d.β), 1), (abs(x - μ) / d.α)^d.β) * inv(2)
        return typeof(v)(1/2) + sign(x - μ) * v
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
    b = 2.0 * rand(Bernoulli()) -1
    return d.μ + inv(sqrt(d.α)) * rand(rng, g)^inv(d.β) * b
end