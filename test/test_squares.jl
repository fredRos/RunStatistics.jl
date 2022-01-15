# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test

@testset "squares" begin
    T_obs = 3.4
    N = 30
    
    @test squares_pvalue(T_obs, N) ≈ 0.818208032939994
    @test squares_cdf(T_obs, N) ≈ 0.18179196706000597
end