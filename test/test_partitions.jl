# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).#

using RunStatistics
using Test

@testset "partitions" begin

    p = RunStatistics.partition(6, 3)

    for i = 0:2

        sum = 0

        for j = 2:(p.h + 1)
               
            sum += p.c[j] * p.y[j]
        end

        @test sum == 6

        RunStatistics.iterate!(p)
    end 
    
    @test RunStatistics.final_partition(p)

end


# NOTE:
#= 

the number of partitions of an integer n into k parts satisfies the recurrence:
p_k(n) = p_k(n - k) + p_k - 1(n - 1)
with p_0(0) = 1 and p_k(n) = 0 uf n <= 0 or k <= 0 and n,k are not both 0

Patition(6, 3)
should therefore yield:

2*1 + 1*4
1*1 + 1*2 + 1*3
3*2

as the partitions (note that the actual output of the function first has to be interpreted correctly, see docstring)

=#
