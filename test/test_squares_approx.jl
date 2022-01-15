# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test

@testset "squares_approx" begin
    
    T_obs = 20
    N = 30
    n = 100
    epsrel = 0.1
    epsabs = 0.001

    @test squares_pvalue_approx(T_obs, N, n, epsrel, epsabs) ≈ 0.24164150728705625
    @test squares_cdf_approx(T_obs, N, n, epsrel, epsabs) ≈ 0.7583584927129438
end