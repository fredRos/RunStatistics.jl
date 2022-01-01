# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).


using RunStatistics
using Test

@testset "squares" begin
    @test pvalue(3.4, 30) == 0.818208032939994
    @test cumulative(3.4, 30) == 0.18179196706000597
end

# test agains the c++ version 