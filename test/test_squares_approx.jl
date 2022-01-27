# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test

@testset "squares_approx" begin
    
    T_obs = 20
    L = 800
    N = 80
    n = 10
    epsp_1 = 10^(-2)
    #TODO: check for zero case for "epsp". throws error in Delta function.
    epsp_2 = 0
       

    @test RunStatistics.squares_pvalue_approx(T_obs, N, n, epsp_1) ≈ 0.07082230261169509
    @test RunStatistics.squares_cdf_approx(T_obs, N, n, epsp_1) ≈ 0.9291776973883049

    @test RunStatistics.squares_cdf_approx(T_obs, L, epsp_1) ≈ 0.9291776973883049
    @test RunStatistics.squares_pvalue_approx(T_obs, L, epsp_1) ≈ 0.07082230261169509


end