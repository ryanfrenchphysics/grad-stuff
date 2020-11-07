module Domain
include("./Lattice.jl")
import .Lattice

lat = Lattice

#=
Define methods and variables related to domain
=#

nx = 0
ny = 0
c = 1.0
ν = 1.0
dim_lu = 0.0
qdim = 9
τ = 0.0
is_interiorsolid = true
is_solid = true

gforcex = 0.0
gforcey = 0.0003

solid = zeros(Int64, (nx, ny))

# Macroscopic velocity
u = zeros(Float64, (nx, ny, 2))
density = zeros(Float64, (nx, ny))
𝑓ᵉ = zeros(Float64, (nx, ny, qdim))   # Distribution function at equilibrium
𝑓ₙ = zeros(Float64, (nx, ny, qdim))   # Distribution function at iteration = n
𝑓ₜ = zeros(Float64, (nx, ny, qdim))   # Distribution function temporary


function getvelocity!(uin, 𝑓in)
    uin = fill!(uin, 0.0)
    for i = 1:nx, j in 1:ny, k in 1:dimq
        uin[i, j, :] .+= 𝑓in[i, j, k] .* lat.directions[k, :] .* c
    end
    uin ./= density
end

#=
function update_macros!(ρin, uin)
    # Clear u and ρ arrays
    fill!(uin, 0.0)
    fill!(ρin, 0.0)

    for i = 1:nx, j in 1:ny
        if is_solidnode
            for a = 1:qdim
                ρ[i, j] .+= 𝑓[i, j, a]
                uin[i, j, 1] .+= ex * 𝑓[i, j, a]
                uin[i, j, 2] .+= ey * 𝑓[i, j, a]
            end
            uin[i, j, 1] ./= ρin[i, j]
            uin[i, j, 2] ./= ρin[i, j]
        end
    end
end
=#


function stream!(𝑓in)
    fill!(𝑓in, 0.0)
    fout = zeros(Float64, (nx, ny, dimq))
    for i in 1:nx, j in 1:ny, k in 1:dimq
        ixn = i + lat.directions[k][1]

        if ixn < 0
            ixn = nx - 1
        elseif ixn >= nx
            ixn = 0
        end

        iyn = j + lat.directions[k][2]

        if iyn < 0
            iyn = ny - 1
        elseif iyn >= ny
            iyn = 0
        end

        𝑓in[ixn, jxn, k] = 𝑓out[i, j, k]
    end
    return 𝑓in
end

function 𝑓eq(ρin, uin)

end


end # Module Domain
