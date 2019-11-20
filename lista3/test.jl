include("module.jl")

using Test

@testset "Bisect tests" begin
    @test bisect(x -> x^2 - 1, -10.0, 10.0, 0.01, 0.1)[4] == 1
    @test bisect(x -> x^2 - 1, 0.0, 2.0, 0.01, 0.1)[1] ≈ 1.0 atol=0.1
end

@testset "Newton tests" begin
end

@testset "Secant tests" begin
end