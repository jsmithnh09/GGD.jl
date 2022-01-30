module GGD

using Statistics, Distributions, SpecialFunctions
import SpecialFunctions: gamma
import Distributions: Gamma, sampler, rand, Bernoulli

export
    ### methods
#    ggdrnd,     # sample from GGD distribution (utilizes Distributions package)
#    ggdmle,     # maximum likelihood function of shape β
#    ggdgcm,     # global convergence method function for β-estimation
#    newton,     # Newton-Raphson estimation routine for root-finding
    params,     # get the distribution parameters (GGD Type)
    shape,      # get the shape parameter (GGD Type)
    scale,      # get the scale parameter (GGD Type)
    location,   # get the value of mu (GGD Type)
    std,        # get the standard deviation (GGD Type)
    mean,       # get the location parameter (GGD Type)
    median,     # median of (GGD Type)
    mode,       # get the mode of the distribution.
    var,        # get the variance of the distribution.
    skewness,   # skewness of the distribution.
    kurtosis,   # EXCESS KURTOSIS, includes the subtraction of 3.
    entropy,    # entropy function of the distribution.
    logpdf,     # log of the PDF function.
    cdf,        # CDF of GGD curve. Uses Gamma Distribution for lower incomplete gamma use.
    pdf,        # PDF of GGD curve.
    rand,       # sample from GGD distribution.



    ### types
    GeneralizedGaussian, # generic type containing params (μ, σ, α, β)

    ### source files

    # types constructor
    include("generalizedgaussian.jl")

    
    """
    A Julia Package for the Generalized Gaussian Distribution (GGD).

    API:
    - `d = GGD(params...)` creates the GGD Type, containing the relevant
        distribution parameters, including mu, sigma, alpha, and beta.
    - `f = GGDFunction(params...)` creates a GGDFunction type, specific to MLE,
        GCM, and Moment-based functions and their associated derivatives.
    - `Z = ggdrnd(n, d)` will generate "n" samples from the specified GGD
        distribution "d".
    - `pdf(d, x)` computes the probability density function at x for GGD
        distribution "d".
    - `cdf(d, x)` computes the cumulative distribution function at x.

    """

    GGD

end # module
