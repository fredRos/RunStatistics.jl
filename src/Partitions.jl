# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

# Implements the algorithm from squares.cxx in https://github.com/fredRos/runs 
# to obtain the partition of n into k parts 
# returns a triple c, y, h with n = \sum_{i=1}^h c_i * y_i

# the vectors for the multiplicities and parts have a buffer value at index 1

mutable struct Partition 


    # some stuff with storage adresses 

    c::Vector{Int}
    y::Vector{Int}
    n::Int
    h::Int
    done::Bool
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



# IMPORTANT: when reading partitions of n into k parts by this algorithm, ignore the first element c[1], y[1] and do not read beyond c[h + 1], y[h + 1]

function Partition(n::Int, k::Int)

    # @argcheck 0 <= k < n

    # multiplicities
    c = zeros(k+1)
    # parts 
    y = zeros(k+1) 
    # number of distinct parts
    h = 0

    y[1] = -1
    c[1] = 1

    #largest part of any partition of n into k parts
    maxPart = n - k +1


    # define intial partition
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



function iterate!(p::Partition)

    # first check if the final partition has already been reached 
    if final_partition(p)
        p.done = true
    else  
        c = p.c
        y = p.y
        h = p.h

        # running index
        i = h
        # calculated part
        k = c[h + 1]
        #calculated remainder
        r = c[h + 1] * y[h + 1]

        # Update step inside the while loop of algorithm Z from
        # A. Zoghbi: Algorithms for generating integer partitions, Ottawa (1993)

        # calculate remainder
        while (y[h + 1] - y[i]) < 2
            k += c[i]
            r += c[i] * y[i]
            i -= 1
        end

        # update current part when it equals 1
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

        # calculate next parts based on remainder left from previous update
        c[i] = k
        r -= c[i] * y[i]
        h = i 

        # update last modified part if it's the remainder
        if (r == y[i])
            c[i] += 1
            h = i - 1
        # add new part with multiplicity 1
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


function final_partition(p::Partition)
    return p.y[p.h + 1] - p.y[2] <= 1 
end 


#=
function inc!(p::Partition)
    if final_partition(p)
        done = true
    else
        println("he")
    end

    return p
end
=#

#find done??













#=

Not really necessary. 
Also, doesn't work, indexing needs to be fixed for Julia 

# Implementation of Algorithm Z from 
#A. Zoghbi: Algorithms for generating integer partitions, Ottawa (1993)

function Partitions(n::Int, k::Int)
    c = zeros(n)
    y = zeros(n)
    parts = zeros(n)
    y[1] = -1
    c[1] = 1
    y[2] = n
    c[2] = 1
    h = 2

    while c[2] != n
        for i in range(2, h)
            parts[i] = c[i] * y[i]        
        end

        i = h - 1
        k = c[h]
        r = c[h] * y[h]

        while (y[h] - y[i]) < 2 
            k += c[i]
            r += c[i] * y[i]
            i -= 1 

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
                i+= 1 
                y[i] = y[i-1] +  1
            end

            c[i] = k
            r -= c[i] * y[i]
            h = i + 1 

            if r == y[i] 
                c[i] += 1 
                h = i
            else 
                y[h] = r
                c[h] = 1
            end    
        end
    end

    return parts
end

=#
