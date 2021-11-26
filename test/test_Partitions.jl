# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).#

using RunStatistics
using Test

@testset "Partitions" begin
    @test RunStatistics.Partitions(5)
end