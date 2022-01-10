# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).


using RunStatistics
using Test

Tobs = 3.4
N = 30

@testset "squares" begin
    @test pvalue(Tobs, N) ≈ 0.818208032939994
    @test cumulative(Tobs, N) ≈ 0.18179196706000597
end

# test agains the c++ version 