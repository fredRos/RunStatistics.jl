# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

"""
    Partition(n::Int, k::Int, h::Int, c::Vector{Int}, y::Vector{Int})

Express the integer partition of `n` into `k` parts in the mulitplicity representation with `n = \\sum_{i=2}^(h + 1) c_i * y_i`. 

(see https://en.wikipedia.org/wiki/Partition_(number_theory))

`h` is the number of distinct parts, `y` an array containing the distinct parts and `c` an array containing their multiplicities.

NOTE: due to the computation of subsequent paritions with the algorithm used in `next_partition!()` the arrays `y` and `c` only hold relevant 
values for the indices [2, h + 1] 

When reading a partition: ignore the first element of `c` and `y` and do not read beyond `c[h + 1]`, `y[h + 1]`!
"""
mutable struct Partition
    n::Int
    k::Int
    h::Int
    c::Vector{Int}
    y::Vector{Int}
end

"""
    init_partition(n::Int, k::Int)

Initiate the first partition of an integer `n` into `k` parts; arguments must satisfy `0 < k <= n`. 
Returns an object of type `Partition`.

The elements in `y[1]` and `c[1]`, of the arrays `y` and `c` containing the distinct parts and their mulitplicities, are buffer values needed for the computation of the next partition 
in `next_partition!()`.

When reading a partition: ignore the first element of `c` and `y` and do not read beyond `c[h + 1]`, `y[h + 1]`: `n = \\sum_{i=2}^(h + 1) c_i * y_i`.
"""
function init_partition(n::Integer, k::Integer)

    @argcheck 0 < k <= n

    c = zeros(k + 1)
    y = zeros(k + 1)
    h = 0

    y[1] = -1
    c[1] = 1

    maxPart = n - k + 1

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

    return Partition(n, k, h, c, y)
end

"""
    next_partition!(p::Partition)

Compute the next partition of `p`, using a modified version of Algorithm Z from *A. Zoghbi: Algorithms for generating integer partitions, Ottawa (1993)*,
https://www.ruor.uottawa.ca/handle/10393/6506. 

The partition `p` is updated in place, saving memory. 
Returns a `boolean` corresponding to whether the final partition has been reached. 
"""
function next_partition!(p::Partition)

    if p.y[p.h+1] - p.y[2] <= 1
        return true
    else
        c = p.c
        y = p.y
        h = p.h

        i = h
        k = c[h+1]
        r = c[h+1] * y[h+1]

        while (y[h+1] - y[i]) < 2
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
            y[i] = y[i-1] + 1
        end

        c[i] = k
        r -= c[i] * y[i]
        h = i

        if (r == y[i])
            c[i] += 1
            h = i - 1

        else
            y[h+1] = r
            c[h+1] = 1
        end

        p.c = c
        p.y = y
        p.h = h
        return false
    end
end

