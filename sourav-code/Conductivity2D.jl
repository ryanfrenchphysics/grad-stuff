using LinearAlgebra
using NLsolve
#using .JuliaNumpyGradient
using PyCall
using PyPlot

np = pyimport("numpy")

hbar = 1.0
m_e = 1.0
kb = 1.0

function dispersion(μ::Float64, kx, ky, a::Float64, t1::Float64, t2::Float64, t3::Float64)

        return(map(abs, μ .+
                (0.5 * t1 * (map(cos, (a * kx)) + map(cos, (a * ky)))) .+
                (t2 * map(cos, (a * kx)) .* map(cos, (a * ky))) .+
                (0.5 * t3 * (map(cos, (2.0 * a * kx)) + map(cos, (2.0 * a * ky))))))

end


function SC_dispersion(Δ::Float64, Ek::Float64)::Tuple{Float64, Float64}
        ret = sqrt(Δ^2 + Ek^2)
        return ret, -ret
end


function quad_dispersion(kx::Float64, ky::Float64)::Float64
        return 0.5 * (kx^2 + ky^2)
end


function SC_dispersion_comp_OAAT(Δ::Float64,  Ek::Float64)::Tuple{Array{Float64,1},Array{Float64,2}}
        M = [
        Ek -Δ;
        -Δ -Ek
        ]
        return eigvals(M), eigvecs(M)
end



function SC_dispersion_comp(Δ::Float64, Ek::Array{Float64,2})::Array{Float64, 2}
        numk = length(Ek[1, :])
        ret = Array{Float64, 2}(undef, numk, numk)
        for i in 1:numk
                for j in 1:numk
                        evals, evecs = SC_dispersion_comp_OAAT(Δ, Ek[i,j])
                        ret[i, j] = maximum(evals)
                end
        end
        return ret
end


function SC_dispersion_comp_new(Δ::Float64, Ek::Array{Float64,1})::Array{Float64, 1}

        ret = Array{Float64, 1}(undef, length(Ek))
        for i in 1:length(Ek)
                evals, evecs = SC_dispersion_comp_OAAT(Δ, Ek[i])
                ret[i] = maximum(evals)
        end
        return ret
end


function gaussian(E::Float64, Ek::Float64, f::Float64)
        return(
        (2.0 / f) * sqrt(log(2.0) / π) * exp(-(4.0 * log(2.0) * (E - Ek)^2) / f^2)
        )
end


function swapcols(vec::Array{Float64, 2})::Array{Float64, 2}
        vec[:, 1], vec[:, 2] = vec[:, 2], vec[:, 1]
        return vec
end

function SC_coherence_comp_OAAT(Δ::Float64, E::Float64, kx::Float64, ky::Float64, μ::Float64, a::Float64, t1::Float64, t2::Float64, t3::Float64)::Float64
        if E < Δ
                χ = 0.0
        else
                χ = sqrt(E^2 - Δ^2)
        end


        evalsk, evecsk = SC_dispersion_comp_OAAT(Δ, χ)

        if evalsk[1] < 0
                evecsk = swapcols(evecsk)
        end

        evalskp, evecskp = SC_dispersion_comp_OAAT(Δ, dispersion(μ, kx, ky, a, t1, t2, t3))

        if evalskp[1] < 0
                evecskp = swapcols(evecskp)
        end

        S = [
        1.0  0.0;
        0.0 -1.0
        ]

        return(
        ((transpose(evecskp) * (S * evecsk))[1, 1])^2
        )
end


function SC_coherence_comp(Δ::Float64, E::Float64, kx::Array{Float64, 2}, ky::Array{Float64, 2}, μ::Float64, a::Float64, t1::Float64, t2::Float64, t3::Float64)::Array{Float64, 2}

        numk = length(kx[1, :])
        ret = Array{Float64, 2}(undef, numk, numk)
        for i in 1:numk
                for j in 1:numk
                        ret[i, j] = SC_coherence_comp_OAAT(Δ, E, kx[i, j], ky[i, j], μ, a, t1, t2, t3)
                end
        end
        return ret
end


function coherence(Δ::Float64, E::Float64, kx::Float64, ky::Float64)::Float64
        return 0.5 * (
        1.0 + (sqrt(E^2 - Δ^2) * dispersion(0.0, kx, ky, 1.0, 1.0, 0.0, 0.0) - Δ^2) /
        (E * SC_dispersion(Δ, dispersion(0.0, kx, ky, 1.0, 1.0, 0.0, 0.0))[1])
        )
end


function GenerateHermite(n)
    Hermite=Function[]

    push!(Hermite,x->1);
    push!(Hermite,x->2*x);

    for ni in 3:n
        push!(Hermite,x->2.0*x*Hermite[ni-1](x)-2.0*(ni - 2)*Hermite[ni-2](x))
    end
    # 2*x*(2x) - (2 * 3 * 1) = 4x^2 - 4
    return Hermite
end


function hermval(x::Float64, cs::Array{Float64, 1})
        # Return c0 * H0(x) + ... + cn * Hn(x)
        ret = cs[1]
        herm = GenerateHermite(length(cs))
        for i in 2:length(cs)
                ret += (cs[i] * map(herm[i], x))
        end
        return ret
end


function integrator_delta_herm(integrate::Array{Float64, 2}, dispersion::Array{Float64, 2}, dk::Float64, E::Float64, n::Int64, f::Float64)::Float64

        coeff = Array{Float64, 1}(undef, n+1)

        for i in 1:length(coeff)
                if ((i-1) % 2) == 0
                        coeff[i] = ((-1)^(i-1)) /
                        (factorial((i-1)) * 4.0^(i-1) *
                        sqrt(π))
                else
                        coeff[i] = 0.0
                end
        end


        displen = length(dispersion)
        exparr = map(exp, (-(((E .- dispersion).^2.) ./ (4.0 * f))))

        hval = Array{Float64, 2}(undef, displen, displen)
        hermarray = E .- dispersion
        for i in 1:displen
                for j in 1:displen
                        hval[i, j] = hermval(hermarray[i, j], coeff)
                end
        end

        func = exparr .* hval .* (1.0 / (2.0 * sqrt(π * f)))
        ret = func .* integrate

        return sum(ret) * dk^2

end

function maskeq(arr::Array{Float64, 2}, val::Float64)::Array{Float64, 2}
        len = length(arr[:, 1])
        mask = Array{Float64, 2}(undef, len, len)
        for i in 1:len
                for j in 1:len
                        if arr[i, j] == val
                                mask[i, j] = 1.0
                        else
                                mask[i, j] = 0.0
                        end
                end
        end
        return mask
end

function maskleq(arr::Array{Float64, 2}, val::Float64)::Array{Float64,1}
        return arr .<= val
end


## DOESN'T LOOK LIKE WE'RE USING THIS
# function integrator_delta_grad(integrate::Array{Float64, 2}, dispersion::Array{Float64, 2}, dk::Float64, E::Float64)::Float64
#
#         mask = maskeq(dispersion, E)
#
#         return 1.0
#
# end


function BZ_integrator(func::Array{Float64, 2}, dk::Float64)::Float64
        return sum(func) .* dk.^2
end


function scattering_SC(Δ::Float64, a::Float64, E::Float64, μ::Float64, t1::Float64, t2::Float64, t3::Float64, numk::Int64, n::Int64, f::Float64)::Float64
        kx = collect(range(-π/a, stop=π/a, length=numk))
        ky = collect(range(-π/a, stop=π/a, length=numk))

        Y = [i for i in kx, j in ky]
        X = [j for i in kx, j in ky]
        dk = 2.0 * π / (a * numk)

        disparr = Array{Float64, 2}(undef, numk, numk)
        for i in 1:numk
                for j in 1:numk
                        disparr[j,i] = dispersion(μ, X[j,i], Y[j,i], a, t1, t2, t3)
                end
        end


        τ = integrator_delta_herm(
                SC_coherence_comp(Δ, E, X, Y, μ, a, t1, t2, t3), SC_dispersion_comp_new(Δ,
                        disparr), dk, E, n, f)

        return τ^(-1)
end


function scattering_n(a::Float64, E::Float64, μ::Float64, t1::Float64, t2::Float64, t3::Float64, numk::Int64, kx::Array{Float64,1}, ky::Array{Float64,1}, n::Int64, f::Float64)::Float64
        dk = 2.0 * π / (a * numk)
        ones = fill!(Array{Float64, 2}(undef, numk, numk), 1.0)

        disparr = Array{Float64, 1}(undef, numk)

        for i in numk
                disparr[i] = dispersion(
                        μ, kx[i], ky[i], a, t1, t2, t3)
        end

        τ = integrator_delta_herm(ones, disparr, dk, E, n, f)
        return τ^(-1)
end


function fermi_deriv(E::Array{Float64,2}, T::Float64)::Array{Float64,2}

        return (
        map(exp, (E ./ (kb .* T))) ./
        ((1.0 .+ map(exp, ((E ./ (kb .* T))))).^2)
        )
end


function FullG!(f, gap, T)
        m = 500
        f[1] = log(T)
        for n in 1:m
                nh = n - 0.5
                sq2 = sqrt((T * nh)^2 + gap[1]^2)
                f[1] += ((1.0 / nh) - (T / sq2))
        end
end


function Delta(T)
        if T < 1.0
                return abs(2.0 * π * nlsolve((F,x) -> FullG!(F,x,T), [0.3]).zero[1])
        end
        return 0.0
end


function p(s, var)
        println(s, ":")
        println(var)
        print("\n\n")
end


function thermal_conductivity(a::Float64, numk::Int64, n::Int64, f::Float64, numT::Int64, Tend::Float64, μ::Float64, t1::Float64, t2::Float64, t3::Float64)::Tuple{Array{Float64,1}, Array{Float64,1}, Array{Float64,1}}

        temps = Array{Float64, 1}(undef, numT)
        kxx = Array{Float64, 1}(undef, numT)
        kxy = Array{Float64, 1}(undef, numT)

        dT = Tend / numT
        T = dT

        kx = collect(range(-π/a, stop=π/a, length=numk))
        ky = collect(range(-π/a, stop=π/a, length=numk))
        # p("kx", kx)
        # p("ky", ky)

        Y = [i for i in kx, j in ky]
        X = [j for i in kx, j in ky]
        # p("X", X)
        # p("Y",Y)

        dk = 2.0 * π / (a * numk)
        disp = Array{Float64, 2}(undef, numk, numk)

        # for i in numk
        #         for j in numk
        #                 disp[i,j] = dispersion(μ, X[i,j], J[i,j], a, t1, t2, t3)
        #         end
        # end


        disp = dispersion(μ, X, Y, a, t1, t2, t3)
        dispSCcomp = Array{Float64, 2}(undef, numk, numk)

        disprows = length(disp[:,1])
        dispcols = length(disp[1,:])
        τ = Array{Float64, 2}(undef, disprows, dispcols)

        for i in 1:disprows
                for j in dispcols
                        τ[i,j] = scattering_n(
                        a, disp[i,j], μ, t1, t2, t3, numk, n, f)
                end
        end

        E² = disp.^2
        # p("Esquared", E²)

        vx = np.gradient(disp, dk)[2]
        vy = np.gradient(disp, dk)[1]
        # p("vx", vx)
        # p("vy", vy)

        for t in 1:numT
                # p("temp num", t)
                println("Temp # ", t, " out of ", numT)
                # println("-------------------\n")
                delta = Delta(T)
                dispSC = SC_dispersion_comp(delta, disp)
                # p("delta",delta)
                # p("dispSC",dispSC)

                dispSCrows = length(dispSC[:,1])
                dispSCcols = length(dispSC[1,:])
                τSC = Array{Float64, 2}(undef, dispSCrows, dispSCcols)
                for i in 1:dispSCrows
                        for j in 1:dispSCcols
                                τSC[j,i] = scattering_SC(delta, a, dispSC[j,i], μ, t1, t2, t3, numk, n, f)
                        end
                end
                # p("TauSC", τSC)

                E²SC = dispSC.^2
                # p("EsquaredSC",E²SC)

                vxSC = np.gradient(dispSC, dk)[2]
                vySC = np.gradient(dispSC, dk)[1]

                # p("vxSC", vxSC)
                # p("vySC", vySC)

                fermi = fermi_deriv(disp, T)
                fermiSC = fermi_deriv(dispSC, T)
                # p("fermi",fermi)
                # p("fermiSC",fermiSC)

                temps[t] = T

                kxxn = BZ_integrator(E² .* τ .* fermi .* vx .* vx, dk) ./ (T^2)
                kxyn = BZ_integrator(E² .* τ .* fermi .* vx .* vy, dk) ./ (T^2)
                kxxSC = BZ_integrator(E²SC .* τSC .* fermiSC .* vxSC .* vxSC, dk) ./ (T^2)
                kxySC = BZ_integrator(E²SC .* τSC .* fermiSC .* vxSC .* vySC, dk) ./ (T^2)

                # p("kxxn", kxxn)
                # p("kxyn", kxyn)
                # p("kxxSC", kxxSC)
                # p("kxySC", kxySC)


                kxx[t] = (kxxSC / kxxn)
                kxy[t] = (kxySC / kxxn)

                T += dT
                # p("Next T", T)
        end
        println(temps)
        println(kxx)
        println(kxy)

        plot(temps, kxx)
        plot(temps, kxy)
        #plot.show()
        return temps, kxx, kxy
end



function thermal_conductivity_OAAT(a::Float64, numk::Int64, n::Int64, f::Float64, T::Float64, μ::Float64, t1::Float64, t2::Float64, t3::Float64, Tc::Float64, Ec::Float64)::Tuple{Array{Float64,1}, Array{Float64,1}, Array{Float64,1}}

        kx = collect(range(-π/a, stop=π/a, length=numk))
        ky = collect(range(-π/a, stop=π/a, length=numk))
        # p("kx", kx)
        # p("ky", ky)

        Y = [i for i in kx, j in ky]
        X = [j for i in kx, j in ky]
        # p("X", X)
        # p("Y",Y)

        dk = 2.0 * π / (a * numk)
        disp = dispersion(μ, X, Y, a, t1, t2, t3)
        Emask = disp .<= Ec
        disp = disp[Emask]
        X = X[Emask]
        Y = Y[Emask]


        disprows = length(disp)
        τ = Array{Float64, 1}(undef, disprows)

        for i in 1:disprows
                τ[i] = scattering_n(
                a, disp[i], μ, t1, t2, t3, numk, X, Y, n, f)
        end

        E² = disp.^2
        # p("Esquared", E²)

        vx = np.gradient(disp, dk)[2]
        vy = np.gradient(disp, dk)[1]
        vx = vx[Emask]
        vy = vy[Emask]
        # p("vx", vx)
        # p("vy", vy)
        delta = Delta(T) * Tc

        dispSC = SC_dispersion_comp_new(delta, disp)

        dispSCrows = length(dispSC)
        τSC = Array{Float64, 1}(undef, dispSCrows)
        for i in 1:dispSCrows
                τSC[i] = scattering_SC(delta, a, dispSC[i], μ, t1, t2, t3, numk, X, Y, n, f)
        end

        E²SC = dispSC.^2

        vxSC = np.gradient(dispSC, dk)[2]
        vySC = np.gradient(dispSC, dk)[1]
        vxSC = vxSC[Emask]
        vySC = vySC[Emask]

        fermi = fermi_deriv(disp, T)
        fermiSC = fermi_deriv(dispSC, T)


        kxxn = BZ_integrator(E² .* τ .* fermi .* vx .* vx, dk) ./ (T^2)
        kxyn = BZ_integrator(E² .* τ .* fermi .* vx .* vy, dk) ./ (T^2)
        kxxSC = BZ_integrator(E²SC .* τSC .* fermiSC .* vxSC .* vxSC, dk) ./ (T^2)
        kxySC = BZ_integrator(E²SC .* τSC .* fermiSC .* vxSC .* vySC, dk) ./ (T^2)

        return T, (kxxSC/kxxn), (kxySC/kxxn)
end
