__precompile__()
module LBMTunnel

include("constants.jl")
using Parameters
using PyPlot
using PyCall
# anim = pyimport("matplotlib.animation")

const anim = PyNULL()

function __init__()
    copy!(anim, pyimport("matplotlib.animation"))
end


export tunnel_main, Distributions, Macros, Barrier, stream!, collide!, boundaryconditions!, gen_circle
#======== Distribtion functions w/ initial values ==========#

@with_kw mutable struct Distributions
    zero = four9ths .* (ones(Float64, δy, δx) .- (1.5 * u₀^2))
    N = one9th .* (ones(Float64, δy, δx) .- (1.5 * u₀^2))
    S = one9th .* (ones(Float64, δy, δx) .- (1.5 * u₀^2))
    E = one9th .* (ones(Float64, δy, δx) .+ (3 * u₀ + 4.5 * u₀^2 - 1.5 * u₀^2))
    W = one9th .* (ones(Float64, δy, δx) .+ (3 * u₀ + 4.5 * u₀^2 - 1.5 * u₀^2))
    NE = one36th .* (ones(Float64, δy, δx) .+ (3 * u₀ + 4.5 * u₀^2 - 1.5 * u₀^2))
    SE = one36th .* (ones(Float64, δy, δx) .+ (3 * u₀ + 4.5 * u₀^2 - 1.5 * u₀^2))
    NW = one36th .* (ones(Float64, δy, δx) .+ (3 * u₀ + 4.5 * u₀^2 - 1.5 * u₀^2))
    SW = one36th .* (ones(Float64, δy, δx) .+ (3 * u₀ + 4.5 * u₀^2 - 1.5 * u₀^2))
end

mutable struct Macros
    ρ::Matrix{Float64}
    ux::Matrix{Float64}
    uy::Matrix{Float64}

    # Constructors:
    function Macros(f::Distributions)
        ρ = (f.zero .+ f.N .+ f.S .+ f.E .+ f.W .+ f.NE .+ f.NW .+ f.SE .+ f.SW)
        ux = (f.E .- f.W .+ f.NE .- f.NW .+ f.SE .- f.SW) ./ ρ
        uy = (f.N .- f.S .+ f.NE .+ f.NW .- f.SE .- f.SW) ./ ρ
        return new(ρ, ux, uy)
    end

    Macros(r, x, y) = new(r, x, y)
end


#=============== Define barrier and BCS ================#
struct Barrier
    barrier::Matrix{Bool}
    N::Matrix{Bool}
    S::Matrix{Bool}
    E::Matrix{Bool}
    W::Matrix{Bool}
    NE::Matrix{Bool}
    NW::Matrix{Bool}
    SE::Matrix{Bool}
    SW::Matrix{Bool}
    BC::Int64

    # Constructor:
    function Barrier(bar, bc)
        bn = circshift(bar, (1,0))
        bs = circshift(bar, (-1,0))
        be = circshift(bar, (0,1))
        bw = circshift(bar, (0,-1))
        bne = circshift(bar, (1,1))
        bnw = circshift(bar, (1,-1))
        bse = circshift(bar, (-1,1))
        bsw = circshift(bar, (-1,-1))
        new(bar, bn, bs, be, bw, bne, bnw, bse, bsw, bc)
    end
end


# Create barrier functions below. Returns Barrier object.
function gen_circle(xoff::Int64, yoff::Int64, rad::Int64, BC::Int64=BC_BOUNCEBACK)
    barr = Array{Bool}(undef, δy, δx)
    for i in 1:δx
        for j in 1:δy
            if (i - xoff)^2 + (j - yoff)^2 <= rad^2
                barr[j , i] = true
            else
                barr[j, i] = false
            end
        end
    end
    return Barrier(barr, BC)
end


function gen_rect(xoff::Int64, yoff::Int64, xlen::Int64, ylen::Int64, BC::Int64=BC_BOUNCEBACK)
    barr = Array{Bool}(undef, δy, δx)
    for i in 1:δx
        for j in 1:δy
            if (xoff - xlen/2) <= i <= (xoff + xlen/2) && (yoff - ylen/2) <= j <= (yoff + ylen/2)
                barr[j, i] = true
            else
                barr[j, i] = false
            end
        end
    end
    return Barrier(barr, BC)
end


# Apply bounce-back boundary conditions
function boundaryconditions!(f::Distributions, b::Barrier)
    for i in 1:δx
        for j in 1:δy
            if b.N[j, i] == true
                f.N[j, i] = f.S[j,i]
            end
            if b.NE[j, i] == true
                f.NE[j, i] = f.SW[j,i]
            end
            if b.NW[j, i] == true
                f.NW[j, i] = f.SE[j,i]
            end
            if b.E[j, i] == true
                f.E[j, i] = f.W[j,i]
            end
            if b.W[j, i] == true
                f.W[j, i] = f.E[j,i]
            end
            if b.S[j, i] == true
                f.S[j, i] = f.N[j,i]
            end
            if b.SE[j, i] == true
                f.SE[j, i] = f.NW[j,i]
            end
            if b.SW[j, i] == true
                f.SW[j, i] = f.NE[j,i]
            end
        end
    end

    # f.N[b.N] .= f.S[b.barrier]
    # f.S[b.S] .= f.N[b.barrier]
    # f.E[b.E] .= f.W[b.barrier]
    # f.W[b.W] .= f.E[b.barrier]
    # f.NE[b.NE] .= f.SW[b.barrier]
    # f.NW[b.NW] .= f.SE[b.barrier]
    # f.SE[b.SE] .= f.NW[b.barrier]
    # f.SW[b.SW] .= f.NE[b.barrier]
end


# Stream (i.e., take step)
function stream!(f::Distributions)
    # Shift every array in given direction
    f.N = circshift(f.N, (1, 0))
    f.S = circshift(f.S, (-1, 0))
    f.E = circshift(f.E, (0, 1))
    f.W = circshift(f.W, (0, -1))

    # f.NE = circshift(f.NE, (1, 0))
    # f.NW = circshift(f.NW, (1, 0))
    # f.SE = circshift(f.SE, (-1, 0))
    # f.SW = circshift(f.SW, (-1, 0))
    f.NE = circshift(f.NE, (1, 1))
    f.NW = circshift(f.NW, (1, -1))
    f.SE = circshift(f.SE, (-1, 1))
    f.SW = circshift(f.SW, (-1, -1))
end


function collide!(f::Distributions, mac::Macros)
    # Update macro values:
    mac.ρ = (f.zero .+ f.N .+ f.S .+ f.E .+ f.W .+ f.NE .+ f.NW .+ f.SE .+ f.SW)
    mac.ux = (f.E .- f.W .+ f.NE .- f.NW .+ f.SE .- f.SW) ./ mac.ρ
    mac.uy = (f.N .- f.S .+ f.NE .+ f.NW .- f.SE .- f.SW) ./ mac.ρ

    # Useful variables for collision step:
    ux2 = mac.ux.^2
    uy2 = mac.uy.^2
    u2 = ux2 .+ uy2
    uxuy = mac.ux .* mac.uy
    A = 1.0 .- 1.5 .* u2

    f.zero = (1.0 .- ω) .* f.zero .+ (ω .* four9ths .* mac.ρ .* A)

    f.N = (1.0 .- ω) .* f.N .+ (ω .* one9th .* mac.ρ .* (A .+ 3 .* mac.uy .+ 4.5 .* uy2))
    f.S = (1.0 .- ω) .* f.S .+ (ω .* one9th .* mac.ρ .* (A .- 3 .* mac.uy .+ 4.5 .* uy2))
    f.E = (1.0 .- ω) .* f.E .+ (ω .* one9th .* mac.ρ .* (A .+ 3 .* mac.ux .+ 4.5 .* ux2))
    f.W = (1.0 .- ω) .* f.W .+ (ω .* one9th .* mac.ρ .* (A .- 3 .* mac.ux .+ 4.5 .* ux2))

    f.NE = (1.0 .- ω) .* f.NE .+ (ω .* one36th .* mac.ρ .* (A .+ 3 .* (mac.ux .+ mac.uy) .+ 4.5 .* (u2 .+ 2 .* uxuy)))
    f.NW = (1.0 .- ω) .* f.NW .+ (ω .* one36th .* mac.ρ .* (A .+ 3 .* (-mac.ux .+ mac.uy) .+ 4.5 .* (u2 .- 2 .* uxuy)))
    f.SE = (1.0 .- ω) .* f.SE .+ (ω .* one36th .* mac.ρ .* (A .+ 3 .* (mac.ux .- mac.uy) .+ 4.5 .* (u2 .- 2 .* uxuy)))
    f.SW = (1.0 .- ω) .* f.SW .+ (ω .* one36th .* mac.ρ .* (A .+ 3 .* (-mac.ux .- mac.uy) .+ 4.5 .* (u2 .+ 2 .* uxuy)))

    # Force rightward flow at ends:

    f.E[:,1] .= one9th .* (1 .+ 3 .* u₀ .+ 3.0 .* u₀.^2)
    f.W[:,1] .= one9th .* (1 .- 3 .* u₀ .+ 3.0 .* u₀.^2)
    f.NE[:,1] .= one36th .* (1 .+ 3 .* u₀ .+ 3.0 .* u₀.^2)
    f.SE[:,1] .= one36th .* (1 .+ 3 .* u₀ .+ 3.0 .* u₀.^2)
    f.NW[:,1] .= one36th .* (1 .- 3 .* u₀ .+ 3.0 .* u₀.^2)
    f.SW[:,1] .= one36th .* (1 .- 3 .* u₀ .+ 3.0 .* u₀.^2)
end

function curl(mac::Macros)
    return circshift(mac.uy, (-1, 1)) - circshift(mac.uy, (1, 1)) - circshift(mac.ux, (-1, 0)) + circshift(mac.ux, (1, 0))
end

function tunnelanimate(i, f::Distributions, mac::Macros, bar_img, fluid_img)
    for step in 1:STEPS_ANIMATION
        stream!(f)
        collide!(f, mac)
    end

    fluid_img.set_array(curl(mac))
    return (fluid_img, bar_img)
end

function tunnelplot(f::Distributions, mac::Macros, b::Barrier)
    tunnelfig = plt.figure()
    fluidimg = plt.imshow(curl(mac), origin="lower", cmap=plt.get_cmap("jet"), interpolation="spline16")
    barrierarr = zeros(Int64, δy, δx, 4)

    for i in 1:δx
        for j in 1:δy
            if b.barrier[j, i] == true
                barrierarr[j, i, 4] = 255
            end
        end
    end
    barrierimg = plt.imshow(barrierarr, origin="lower", interpolation="quadric")

    # plt.title(ANIMATION_TITLE)
    plt.subplots_adjust(top=0.88)

    animate = anim.FuncAnimation(tunnelfig, tunnelanimate, fargs=(f, mac, barrierimg, fluidimg), interval=1, blit=true, repeat=true)
    # plt.clim(vmin=0, vmax=1)
    plt.colorbar(fluidimg)
    plt.show()
end

function tunnel_main()
    distros = Distributions()
    macros = Macros(distros)
    # barrier = gen_circle(100,round(Int,δy/2),15)
    barrier = gen_rect(100, round(Int, δy/2), 30, 50)

    tunnelfig = plt.figure()
    fluidimg = plt.imshow(curl(macros), origin="lower", cmap=plt.get_cmap("jet"), interpolation="spline16")
    barrierarr = zeros(Int64, δy, δx, 4)

    for i in 1:δx
        for j in 1:δy
            if barrier.barrier[j, i] == true
                barrierarr[j, i, 4] = 255
            end
        end
    end
    barrierimg = plt.imshow(barrierarr, origin="lower", interpolation="quadric")

    # plt.title(ANIMATION_TITLE)
    plt.subplots_adjust(top=0.88)


    function nextframe(i)
        for step in 1:STEPS_ANIMATION
            stream!(distros)
            boundaryconditions!(distros, barrier)
            collide!(distros, macros)
        end

        fluidimg.set_array(curl(macros))
        return (fluidimg, barrierimg)
    end

    animate = anim[:FuncAnimation](tunnelfig, nextframe, interval=1, blit=true)
    plt.clim(vmin=0, vmax=1)
    plt.colorbar(fluidimg)
    plt.show()


    # tunnelplot(distros, macros, barrier)

end


end # module LBMTunnel
