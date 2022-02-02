# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test

# Before you make a release, you should copy over all the tests I implemented in https://github.com/fredRos/runs/tree/master/test and make sure the numbers agree. There may some discrepancy because a different numerical integration is used but otherwise you should get the same numbers

@testset "squares" begin
    T_obs = 3.4
    N = 30

    @test squares_pvalue(T_obs, N) ≈ 0.818208032939994
    @test squares_cdf(T_obs, N) ≈ 0.18179196706000597
end