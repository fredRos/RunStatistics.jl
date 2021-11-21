# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test

@testset "Tobs" begin
    @test RunStatistics.Tobs([-1, 1, 3, -2], 0, 1) == 10
end
