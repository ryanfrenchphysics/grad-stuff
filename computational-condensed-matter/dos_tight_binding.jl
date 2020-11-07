using PyCall
using TimerOutputs

const to = TimerOutput()

py"""
from itertools import permutations
def permut(x, r):
    return list(permutations(x, r))
"""
plt = pyimport("matplotlib.pyplot")

function dispersion(k_perms::Array{Float64, 1}, dimensions::Int64, num_k::Int64,
    hop::Float64, spacing::Float64, E₀::Float64)
    
    val::Float64 = E₀
    @timeit to "dispersion loop" begin
        for i in 1:dimensions
            val -= 2 * hop * cos(k_perms[i] * spacing)
        end
    end
    return val
end

function k_perms(k_mat, dimensions::Int64, num_k::Int64)
    k_vals = Float64[]

    @timeit to "k_perms loop 1" begin
        for i in 1:num_k
            push!(k_vals, k_mat[i][1])
        end
    end

    @timeit to "k_perms gen permut" begin
        k_permutations = py"permut"(k_vals, dimensions)
        k_permutations = [collect(x) for x in k_permutations]
    end

    @timeit to "k_perms vcat" begin
        for i in 1:num_k
            vcat(k_permutations, k_mat[i])
        end
    end

    return k_permutations
end


function gen_E(start::Float64, final::Float64, n::Int64)
    step::Float64 = (final - start) / n

    Es = zeros(Float64, n+1)

    @timeit to "gen_E loop" begin
        for i in 1:n+1
            Es[i] = start + ((i - 1) * step)
        end
    end

    return Es
end


function gen_D_mat(Earr::Array{Float64, 1}, disparr::Array{Float64, 1}, error::Float64, Eₘₐₓ::Float64)
    D = zeros(length(Earr))
    for i in 1:length(Earr)
        states_at_E = 0
        for j in 1:length(disparr)
            if abs(Earr[i] - disparr[j]) <= error
                states_at_E += 1
            end
        end
        D[i] = (states_at_E * Eₘₐₓ)
    end
    return D
end

function DOS(dimensions::Int64, error::Float64, numk::Int64, numE::Int64,
    hop::Float64, spacing::Float64)

    E₀ = convert(Float64, (2.0 * dimensions))

    kmat = zeros(Float64, numk, dimensions)
    stepsize = ((2 * π) / spacing) / numk

    @timeit to "gen k matrix" begin
        for i in 1:numk
            kmat[i,:] .= (-π / spacing) + ((i - 1) * stepsize)
        end
    end

    kperm_mat = k_perms(kmat, dimensions, numk)

    dispersionarray = zeros(length(kperm_mat))

    for i in 1:length(kperm_mat)
        dispersionval = dispersion(kperm_mat[i], dimensions, numk, hop, spacing, convert(Float64, E₀))
        dispersionarray[i] = dispersionval
    end

    Es = gen_E(0.0, convert(Float64, 4.0 * dimensions), numE)
    Eₘₐₓ = maximum(Es)

    @timeit to "gen_D_mat" begin
        D = gen_D_mat(Es, dispersionarray, error, Eₘₐₓ)
    end

    Dₙₒᵣₘ = zeros(length(D))
    for i in 1:length(D)
        Dₙₒᵣₘ[i] = convert(Float64, (D[i] * ((spacing / π)^dimensions)))
    end

    Eₙₒᵣₘ = zeros(length(Es))
    for i in 1:length(Es)
        Eₙₒᵣₘ[i] = (Es[i] / Eₘₐₓ)
    end

    print(to)

    #print(Eₙₒᵣₘ)
    #print(D)
    plt.plot(Eₙₒᵣₘ, Dₙₒᵣₘ)
    plt.title("DOS: $dimensions Dimensions")
    plt.xlabel("Energy")
    plt.ylabel("n")
    plt.show()
end