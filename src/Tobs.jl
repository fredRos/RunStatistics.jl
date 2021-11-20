# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

function Tobs(X::AbstractArray, μ::Real, σ2::Real)
    χ2 = []
    χi = 0
    for i in range(1, length(X))
        
        X[i] > μ ?  χi += (X[i] - μ)^2 / σ2 : (append!(χ2, χi) ; χi = 0)

        #if X[i] > μ
        #    χi += (X[i] - μ)^2 / σ2
        #else
        #    append!(χ2, χi)
       #     χi = 0
       # end 
    end
 
    return maximum(χ2)

end
 

println(Tobs([-1, 1, 3, -2], 0, 1))