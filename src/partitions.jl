# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).


"""
    Partition(c::Vector{Int}, y::Vector{Int}, n::Int, h::Int, done::Bool)

Represent the integer partition of `n` into `k` parts, with n = \\sum_{i=1}^h c_i * y_i; `h` is the number of distinct parts, `y` the parts and `c` their multiplicities.

When reading: ignore first element of `c` and `y` and do not read beyond `c[h + 1]`, `y[h + 1]`.
"""
mutable struct Partition 
    c::Vector{Int}
    y::Vector{Int}
    n::Int
    h::Int
    done::Bool
end

"""
    partition(n::Int, k::Int)

Initiate the first partition of an integer `n` into `k` parts, arguments must satisfy `0 < k <= n`. 

Returns an object of type `Partition`
"""
function partition(n::Int, k::Int)

    @assert 0 < k <= n

    c = zeros(k+1)
    y = zeros(k+1) 
    h = 0

    y[1] = -1
    c[1] = 1

    maxPart = n - k +1

    if (k == 1 || k == n)
        y[2] = maxPart
        c[2] = n / maxPart
        h = 1
        done = true
    else
        y[2] = 1 
        c[2] = k - 1
        y[3] = maxPart
        c[3] = 1
        h = 2
        done = false 
    end

     
    return Partition(c, y, n, h, done)
end 

"""
    iterate!(p::partition)

Compute the next partition of `p`, updating it in place. Returns an object of type `partition`. 
"""
function iterate!(p::Partition)

    if final_partition(p)
        p.done = true
    else  
        c = p.c
        y = p.y
        h = p.h

        i = h
        k = c[h + 1]
        r = c[h + 1] * y[h + 1]

        while (y[h + 1] - y[i]) < 2
            k += c[i]
            r += c[i] * y[i]
            i -= 1
        end

        if c[i] == 1
            if i != 1
                r += c[i] * y[i]
                y[i] += 1
            else
                i = 2
                y[i] = 1
            end
        else
            c[i] -= 1
            r += y[i]
            i += 1
            y[i] = y[i - 1] + 1 
        end

        c[i] = k
        r -= c[i] * y[i]
        h = i 

        if (r == y[i])
            c[i] += 1
            h = i - 1

        else
            y[h + 1] = r
            c[h+ 1] = 1
        end

        p.c = c
        p.y = y
        p.h = h
        return p 
    end
end

"""
    final_partition(p::partition)

Check whether `p` is the final partition. Returns a `boolean`.
"""
function final_partition(p::Partition)
    return p.y[p.h + 1] - p.y[2] <= 1 
end 






#=
struct AbstractPartitionGenerator

    function bool(p::AbstractPartitionGenerator)
        return !p.done
    end

    # some stuff again with adresses
    AbstractPartitionGenerator(p::Partition)

    function final_partition(p::partition)
        return p.y[p.h + 1] - p.y[2] <= 1 
    end
    
    
    # might be a problem, is 'virtual bool' in c++
    done::Bool
    p::typeof(p) # does this make sense?
end



struct PartitionsGenerator <: AbstractPartitionGenerator
    PartitionsGenerator(n::Int, k::Int)

    operator++()::PartitionsGenerator
     # problem with final_partition() in parent, how does the override work in julia?


end


function PartitionGenerator.final_partition()
    return
end

=#