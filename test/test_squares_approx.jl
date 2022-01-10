# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).


using RunStatistics
using Test

Tobs = 20
N = 30
n = 100
epsrel = 0.1
epsabs = 0.001


@testset "squares_approx" begin
    @test approx_pvalue(Tobs, N, n, epsrel, epsabs) ≈ 0.24164150728705625
    @test approx_cumulative(Tobs, N, n, epsrel, epsabs) ≈ 0.7583584927129438
end

