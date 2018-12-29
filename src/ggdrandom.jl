function ggdrnd1(μ, σ, β)
    # determining GGD scale.
    α = σ * sqrt( gamma((1. / β)) / gamma((3. / β)) )

    # generating the distribution.
    d = Gamma((1. / β), 1)

    # generating sampler object
    s = sampler(d)

    # generating gamma sample.
    γ = rand(s)
    Γ = γ .^ ((1. / β))

    # Bernoulli value.
    B = (rand(1)[1] < 0.5)*2 - 1

    # GGD value. Z = μ + 1/√(α) * |Y|¹⁻ᴮ * B
    Z = μ + ((1. / √(α)) * Γ * B)
    return Z
end

function ggdrnd(n, μ, σ, β)
    # pre-allocating empty Array.
    Z̄ = Array{Float64}(undef, n)

    # add the elements to the vector.
    for i = 1:n
        Z̄[i] = ggdrnd1(μ, σ, β)
    end
    return Z̄
end

Z = ggdrnd(100, 0, 1, 2)
