using TimerOutputs

const to = TimerOutput()

function permute!(a::Array{Float64, 1}, r::Int64=0)
    n::Int64 = length(a)
    done = false
    if r == 0
        r = n
    end
    edge::Int64 = r

    # Find index j in (r+1...n) where aⱼ > a_edge
    j::Int64 = edge + 1
    @timeit to "First edge check" begin
    while j <= n && a[edge] >= a[j]
        j += 1
    end
    end

    if j <= n
        # Swap a_edge, aⱼ
        @timeit to "if block" begin
        a[edge], a[j] = a[j], a[edge]
        end
    else
        @timeit to "else block" begin
            # Reverse aᵣ₊₁ to aₙ
            @timeit to "First reverse" begin
                a[edge+1:n] = a[n:-1:edge+1]

            end

            # Find rightmost ascent to left of edge
            i::Int64 = edge - 1
            @timeit to "Finding right ascent" begin
                while i > 0 && a[i] >= a[i + 1]
                    i -= 1
                end
            end

            if i <= 0
                # No more permutations
                done = true
            end

            # Find j in (n ... i+1) where aⱼ > aᵢ
            j = n
            @timeit to "Find j for i" begin
                while j > i && a[i] >= a[j]
                    j -= 1
                end
            end

            # Swap aᵢ and aⱼ
            @timeit to "Swap btwn reversal" begin
                a[i], a[j] = a[j], a[i]
            end

            # Reverse a_(i + 1) to a_(n)
            @timeit to "Last reversal" begin
                a[i+1:n] = a[n:-1:i+1]
            end
        end
    end
    if done == true
        a = []
    end
end

function permutation(a::Array{Float64, 1}, r::Int64=0)
    n::Int64 = length(a)
    if r == 0
        r = n
    end

    next_permutation = sort(a)
    returnarr = [a[1:r]]
    # print(returnarr)
    # print("\n")

    while true
        try
            permute!(next_permutation, r)
        catch err
            if isa(err, BoundsError)
                break
            end
        end
        # print(next_permutation[1:r])
        # print("\n")
        push!(returnarr, next_permutation[1:r])
    end
    print(to)
end
