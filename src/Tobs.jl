# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

#Think about implementing a special case with μ and σ as vectors, so each obsevation X[i] can have an individual normal distribution

function Tobs(X::AbstractArray, μ::Real, σ2::Real)
    χ2 = []
    χi = 0
    for i in range(1, length(X))
        X[i] > μ ?  χi += (X[i] - μ)^2 / σ2 : (append!(χ2, χi); χi = 0)
    end 
    return maximum(χ2)
end
