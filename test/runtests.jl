# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using Test
import RunStatistics
import Documenter

@testset "Package RunStatistics" begin

    include("test_partitions.jl")
    include("test_tobs.jl")
    include("test_squares.jl")
    include("test_squares_approx.jl")

    # doctests
    Documenter.DocMeta.setdocmeta!(
        RunStatistics,
        :DocTestSetup,
        :(using RunStatistics);
        recursive = true
    )

    Documenter.doctest(RunStatistics)

end # testset


