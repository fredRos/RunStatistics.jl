# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

import Test
import RunStatistics
import Documenter

Test.@testset "Package RunStatistics" begin
    include("test_hello_world.jl")
    include("test_Partitions.jl")
    include("test_Tobs.jl")
    # doctests
    Documenter.DocMeta.setdocmeta!(
        RunStatistics,
        :DocTestSetup,
        :(using RunStatistics);
        recursive=true,
    )
    Documenter.doctest(RunStatistics)
end # testset


