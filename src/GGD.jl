"""
Generalized Gaussian Distribution package.
"""
module GGD

using Distributions
using SpecialFunctions: gamma
import Base: rand
import Statistics: std, var, mean
import Random: GLOBAL_RNG, AbstractRNG

### source files
include("generalizedgaussian.jl")

### methods
export
    GeneralizedGaussian,     # principle distribution type.
    params,                  # returns shape, location, and scale parameters.
    shape,                   # shape parameter β.
    scale,                   # scale parameter α, which is related to variance via the shape parameter and gamma functions.
    location,                # location parameter μ, the same as the mean/median/mode.
    std,                     # standard deviation σ.
    mean,                    # mean of the distribution (equal to μ)
    median,                  # median of distribution
    mode,                    # mode of distribution
    var,                     # variance of distribution, (σ²)
    skewness,                # skewness of Generalized Gaussian Distribution is zero.
    kurtosis,                # Excess kurtosis of distribution, (this subtracts 3.)
    entropy,                 # entropy of distribution
    logpdf,                  # log probability density
    cdf,                     # cummulative distribution function
    pdf,                     # probability density function
    gcmsearch,               # global convergence method for estimating ̂β
    rand                     # random sample from the distribution

"""
A Julia Package for the Generalized Gaussian Distribution (GGD).

API:
- `d = GGD(params...)` creates the GGD Type, containing the relevant
    distribution parameters, including mu, sigma, alpha, and beta.
- `z = rand(rng, d)` will generate a sample from the specified GGD
    distribution "d".
- `pdf(d, x)` computes the probability density function at x for GGD
    distribution "d".
- `cdf(d, x)` computes the cumulative distribution function at x.

"""


end # module
