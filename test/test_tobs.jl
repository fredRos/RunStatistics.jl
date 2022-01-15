# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test

@testset "t_obs" begin
    @test RunStatistics.t_obs([-1, 1, 3, -2], 0, 1) == 10
end