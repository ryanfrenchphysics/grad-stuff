module Ex
    using PyPlot
    using LinearAlgebra

    export fromfunction
    export equilibrium
    # from LBM_ex.py in ~/


    #Flow Definition
    const maxiter = 1000    # iterations
    const Re = 220.0        # Reynold's number

    const nx = 520          #
    const ny = 180          # Lattice dimensions/
    const ly = ny - 1.0     # populations
    const q = 9             #

    const cx = nx / 4       #
    const cy = ny / 2       # Cylinder coordinates
    const r = ny / 2        #

    const u_lb = 0.04
    const nu_lb = u_lb * r / Re
    const ω = 1.0 / (3.0 * nu_lb + 0.5) # Relaxation parameter

    # Lattice constants

    # Lattice velocities
    c = [(x, y) for x in [0.0, -1.0, 1.0] for y in [0.0, -1.0, 1.0]]
    t = 1.0 / 36.0 * ones(q)

    noslip = [findall(x->x==.-val, c) for val in c]

    # Unknowns on
    i1 = Int64[]    # Right wall
    i2 = Int64[]    # Vertical middle
    i3 = Int64[]    # Left wall

    # Lattice weights
    for (idx, val) in enumerate(c)
        if norm(val) < 1.1
            t[idx] = 1.0 / 9.0
        end
        if val[1] < 0
            push!(i1, idx)
        elseif val[1] == 0
            push!(i2, idx)
        elseif val[1] > 0
            push!(i3, idx)
        end

    end
    t[1] = 4.0 / 9.0

    function fromfunction(func::Function, shape, type)
        arr = Array{type}(undef, shape)
        for I in CartesianIndices(arr)
            arr[I] = func(Tuple(I)...)
        end
        return arr
    end

    function equilibrium(ρ, u)
        # u is permutims of vel
        u_trans = permutedims(u, (3,2,1)) #(2, 180, 520)
        cu = Array{Float64}(undef, q, size(u_trans)[2], size(u_trans)[3])
        for i in 1:size(u)[1]
            cu[:,:,i] = c * u_trans[:,:,i]
        end
        cu = permutedims(cu, (3,2,1))
        cu = 3.0 .* cu
        usqr = (3.0 / 2.0) * (u[:,:,1].^2 + u[:,:,2].^2)
        feq = zeros(nx, ny, q)

        for i in 1:q
            feq[:,:,i] = ρ .* t[i] .* (1.0 .+ cu[:,:,i] .+ 0.5 .* cu[:,:,i].^2 .- usqr)
        end
        return feq
    end

    vel = vel = fromfunction((x, y,d) -> (1-(d-1))*u_lb*(1.0+1e-4*sin((y-1)/ly*2*pi)),(nx,ny,2), Float64)

end
