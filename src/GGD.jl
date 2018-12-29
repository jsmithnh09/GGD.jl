module GGD

using Statistics, Distributions, SpecialFunctions

import SpecialFunctions: gamma
import Base: rand, randn
import Distributions: Gamma, sampler, rand

export
    ### methods
    ggdrnd,     # sample from GGD distribution (utilizes Distributions package)
    ggdcdf,     # CDF of GGD curve
    ggdpdf,     # PDF of GGD curve
    ggdmle,     # maximum likelihood function of shape β
    ggdgcm,     # global convergence method function for β-estimation
    newton,     # Newton-Raphson estimation routine for root-finding
    params,     # get the distribution parameters (GGD Type)
    beta,       # get the shape parameter (GGD Type)
    alpha,      # get the scale parameter (GGD Type)
    mean,       # get the location parameter (GGD Type)

    ### types
    GGD,         # generic type containing params (μ, σ, α, β)
    GGDFunction # Function type (MLE, GCM, MoLC, Mallat's Method, etc.)

    # TODO: Check Distributions/src/univariate/inversegaussian.jl for proper
    # constructor overview.

    ### source files

    # random sample generator
    include("ggdrandom.jl")

    # estimation functions for β
    include("ggdest.jl")
    include("ggdroot.jl")

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
