# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test

@testset "hello_world" begin
    @test RunStatistics.hello_world() == 42
end
