using GGD
using Test

@testset "Default Distribution" begin
    # default args should generate a default Normal distribution N(0, 1).
    d = GeneralizedGaussian()
    @test mean(d) ≈ 0
    @test std(d) ≈ 1
    @test var(d) ≈ 1
    @test shape(d) ≈ 2
end

d = GeneralizedGaussian(1.8) # going to test with shape = 1.8.
@testset "GGD sampler" begin
    @test_nowarn X = rand(d) # singular case 
    @test_nowarn X = rand(d, (100, 1)) # matrix case
    @test_nowarn X = rand(d, 100) # vector case
end

X = rand(d, 7_000) # higher N = closer shape estimate.
@testset "GCM Search" begin
    @test_nowarn β = gcmsearch(X, 1.0, 100)
    @test_nowarn β = gcmsearch(X) # no args should use absolute moment initial estimate
end

@testset "GCM Interval" begin
    β = gcmsearch(X)
    ci = gcmci(β, length(X))
    @test ci[1] ≤ β ≤ ci[2]
    @test isapprox(shape(d), β, rtol=(ci[2]-ci[1])/ci[2])
end




