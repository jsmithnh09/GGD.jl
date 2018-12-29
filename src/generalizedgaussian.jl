"""
    GeneralizedGaussian( μ, σ, β )

The "Generalized Normal" or "Generalized Gaussian" with shape parameter 'β',
scale 'α', and location 'μ' has probability density function

```math
f(x; \\mu \\alpha \\beta ) = \\frac{\\beta}{2\\alpha\\Gamma(1/\\beta)}
\\text{exp}\\left( -|x - \\mu|/\\alpha \\right),

\\alpha = \\sigma \\sqrt{ \\frac{\\Gamma(1/\\beta)}{\\Gamma(3/\\beta)} }

```

```julia
GeneralizedGaussian(m, s, b)    # Generalized Gaussian with mu, sigma, and beta.
params(d)                       # Get the parameters (m,s,b)
location(d)                     # returns location mu
std(d)                          # returns sigma
scale(d)                        # returns the scale from the sigma/beta function
shape(d)                        # returns the shape beta
```

External Links
[Generalized Normal Distribution on Wikipedia](https://en.wikipedia.org/wiki/Generalized_normal_distribution)

References
[1] Nadarajah, Saralees (September 2005). "A generalized normal distribution".
    Journal of Applied Statistics. 32 (7): 685&ndash, 694.
"""
struct GeneralizedGaussian{T<:Real}
    μ::T
    σ::T
    β::T
    α::T

    function GeneralizedGaussian{T}(μ::T, σ::T, β::T) where T
        σ > zero(σ) || error("standard deviation must be positive.")
        β > zero(β) || error("β must be greater than zero.")
        new{T}(μ, σ, β)
    end
end

GeneralizedGaussian(μ::T, σ::T, β::T) where {T<:Real} = GeneralizedGaussian{T}(μ, σ, β)
GeneralizedGaussian(μ::Real, σ::Real, β::Real) = GeneralizedGaussian(promote(μ, σ, β)...)

function GeneralizedGaussian(μ::Integer, σ::Integer, β::Integer)
    GeneralizedGaussian(Float64(μ), Float64(σ), Float64(β))
end

### Conversion
function convert(::Type{GeneralizedGaussian{T}}, μ::Real, σ::Real, β::Real) where T <: Real
    GeneralizedGaussian(T(μ), T(σ), T(β))
end
function convert(::Type{GeneralizedGaussian{T}}, d::GeneralizedGaussian{S}) where {T <: Real, S <: Real}
    GeneralizedGaussian(T(d.μ), T(d.σ), T(d.β))
end

### Parameter
shape(d::GeneralizedGaussian) = d.β
location(d::GeneralizedGaussian) = d.μ
std(d::GeneralizedGaussian) = d.σ
scale(d::GeneralizedGaussian) = d.σ * sqrt( gamma(1. / d.β) / gamma(3. / d.β) )
params(d::GeneralizedGaussian) = (d.μ, d.σ, d.β)
@inline partype(d::GeneralizedGaussian{T}) where {T<:Real} = T


### Statistic Functions
function mean(d::GeneralizedGaussian{T}) where T<:Real
    (μ, σ, β) = params(d)
    return μ
end

function median(d::GeneralizedGaussian)
    (μ, σ, β) = params(d)
    return μ
end

function mode(d::GeneralizedGaussian)
    (μ, σ, β) = params(d)
    return μ
end

function var(d::GeneralizedGaussian{T}) where T<:Real
    (μ, σ, β) = params(d)
    return σ * σ
end

@inline function skewness(d::GeneralizedGaussian{T}) where T<:Real return 0 end

function kurtosis(d::GeneralizedGaussian{T}) where T<:Real
    (μ, σ, β) = params(d)
    return ( ( gamma(5. / β) * gamma(1. / β) ) / ( gamma(3. / β)^2 ) ) - 3
end

function entropy(d::GeneralizedGaussian{T}) where T<:Real
    (μ, σ, β) = params(d)
    α = scale(d)
    return (1. / β) - log( β / (2 * α * gamma(1. / β)))
end


### Evaluation
function pdf(d::GeneralizedGaussian{T}, x::Real) where T<:Real
    if x == -Inf || x == Inf
        return zero(T)
    else
        (μ, σ, β) = params(d)
        α = scale(d)

        A = β / (2 * α * gamma(1. / β) )
        return A * exp( -(abs(x - μ) / α)^β )
    end
end

function cdf(d::GeneralizedGaussian{T}, x::Real) where T<:Real
    if x == -Inf || x == Inf
        return zero(T)
    else
        (μ, σ, β) = params(d)
        α = scale(d)
        #TODO: construct the lower incomplete gamma function.
        Q = 0.5
        return 0.5 + sign(x - μ)*( Q / (2*gamma(1. / β)) )
    end
end
