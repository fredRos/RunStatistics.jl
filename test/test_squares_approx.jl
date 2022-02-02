# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test

@testset "squares_approx" begin

    T_obs = 20
    N = 30
    n = 100
    epsp = 10^(-2)

    @test squares_pvalue_approx(T_obs, N, n, epsp) ≈ 0.24164150728705625
    @test squares_cdf_approx(T_obs, N, n, epsp) ≈ 0.7583584927129438
end

# I measured the timing for N=96 with the exact and approximate solution. See if julia is really faster. It seems some experiment with threading was done but without threads, I'd be surprised if it's faster than the C++ code built in release mode