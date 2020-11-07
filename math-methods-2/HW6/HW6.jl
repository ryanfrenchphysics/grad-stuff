##############################################################################
# PHSX 567 homework assignment 6: Advection of Atmospheric Aerosols
#
#
# Developed and confirmed execution on:
# OS:                 Manjaro Linux 18.1.5
# Kernel:             Linux 4..19.96-1-MANJARO
# Architecture:       x86-64
# Julia Relesase:     1.5.0-DEV
#
#
#
# Created by: Ryan French
# Date: 04-22-2020 #############################################################################

using DifferentialEquations
using ODEInterfaceDiffEq    # For radau() method
using PyPlot
using ArgParse


#======== CONSTANTS ====================#
const ρ = 1.0                 # Density of air
const ν = 1000.0              # Density of droplet
const w₀ = 10.0               # Speed of vortex
const vortex_rad = 1000.0     # Radius of vortices (m)
const η = 2e-5                # Air viscosity
const g = 10.0                # Gravitational acceleration

const x0 = 1e-30  # 0.0 breaks throws an error...I believe this is a module bug, I'm submitting a report
const xdot0 = 0.0
const y0 = 1500.0
const ydot0 = 0.0

const τ = (2 * π * vortex_rad) / (w₀)  # turnover time
const TMIN = 0.0
const TMAX = 10.0 * τ

const diameter_vals = [
    1e-6, 1e-5, 1e-4, 1e-3
]


s = ArgParseSettings()
@add_arg_table s begin
    "--noplot", "-n"
        help = "Stops plots from showing"
        action = :store_true
    "--scaled-plots", "-s"
        help = "Turns on shared-y axis plots"
        action = :store_true
    "--choose", "-c"
        help = "Choose to run only one part of the homework assignment. Accepts values 1, 2. Default runs both parts."
        arg_type = Int
        default = 0
end


"""
    Take the eddy differential equation and break it into 4
    first-order differential equations
"""
function vectorfield!(du, u, p, t)
    x, vx, y, vy = u

    d = p[1]
    ww = (π * w₀) / 2.0
    ω = π / (2.0 * vortex_rad)
    wx = -ww * cos(ω * x) * cos(ω * y)
    wy = -ww * sin(ω * x) * sin(ω * y)

    ux = wx - vx
    uy = wy - vy

    reynoldx = (ρ * d / η) * abs(ux)
    reynoldy = (ρ * d / η) * abs(uy)
    CDx = (24 / reynoldx) * (1 + ((1 / 6) * reynoldx^(2/3)))
    CDy = (24 / reynoldy) * (1 + ((1 / 6) * reynoldy^(2/3)))
    A = π * (d / 2)^2
    V = (4 / 3) * π * (d / 2)^3
    m = ν * V

    du[1] = vx
    du[2] = (1 / (2 * m)) * CDx * A * ρ * abs(ux) * ux
    du[3] = vy
    du[4] = -g + (1 / (2 * m)) * CDy * A * ρ * abs(uy) * uy
end

"""
    Create subplots for each
"""
function plot4(s, d, type="time", noplot=false, share_y=false)
    if noplot == false
        plt.clf()
        fig, ax = subplots(2, 2, sharex=true, sharey=share_y)
        fig.suptitle("Position of Droplet in Vortex", fontsize=20)

        len = length(s)
        ts = []
        xs = []
        ys = []
        for i in 1:len
            push!(ts, s[i].t)
            xtmp = []
            ytmp = []
            for j in 1:length(s[i].t)
                push!(xtmp, s[i].u[j][1])
                push!(ytmp, s[i].u[j][3])
            end
            push!(xs, xtmp)
            push!(ys, ytmp)
        end

        if type == "time"
            ax[1, 1].plot(ts[1], xs[1], "r", label="x(t)")
            ax[1, 1].plot(ts[1], ys[1], "b", label="y(t)")
            ax[1, 1].set_title("Position of d=$(d[1])", fontsize=16)

            ax[2, 1].plot(ts[2], xs[2], "r", label="x(t)")
            ax[2, 1].plot(ts[2], ys[2], "b", label="y(t)")
            ax[2, 1].set_title("Position of d=$(d[2])", fontsize=16)

            ax[1, 2].plot(ts[3], xs[3], "r", label="x(t)")
            ax[1, 2].plot(ts[3], ys[3], "b", label="y(t)")
            ax[1, 2].set_title("Position of d=$(d[3])", fontsize=16)

            ax[2, 2].plot(ts[4], xs[4], "r", label="x(t)")
            ax[2, 2].plot(ts[4], ys[4], "b", label="y(t)")
            ax[2, 2].set_title("Position of d=$(d[4])", fontsize=16)

            handles, labels = plt.gca().get_legend_handles_labels()
            fig.legend(handles, labels, loc="center right", fancybox=true, fontsize=12, markerscale=2.0)

            fig.text(0.5, 0.04, "time (s)", ha="center", va="center", fontsize=16)
            fig.text(0.06, 0.5, "position (m)", ha="center", va="center", rotation="vertical", fontsize=16)

        elseif type == "position"
            ax[1, 1].plot(xs[1], ys[1], "r")
            ax[1, 1].set_title("Position of d=$(d[1])", fontsize=16)

            ax[2, 1].plot(xs[2], ys[2], "b")
            ax[2, 1].set_title("Position of d=$(d[2])", fontsize=16)

            ax[1, 2].plot(xs[3], ys[3], "y")
            ax[1, 2].set_title("Position of d=$(d[3])", fontsize=16)

            ax[2, 2].plot(xs[4], ys[4], "g")
            ax[2, 2].set_title("Position of d=$(d[4])", fontsize=16)

            # handles, labels = plt.gca().get_legend_handles_labels()
            # fig.legend(handles, labels, loc="center right", fancybox=true, fontsize=12, markerscale=2.0)

            fig.text(0.5, 0.04, "horizontal position (m)", ha="center", va="center", fontsize=16)
            fig.text(0.06, 0.5, "vertical position (m)", ha="center", va="center", rotation="vertical", fontsize=16)
        end
        plt.grid(false)
        plt.show()
    end
end


function optimsize(initlow, inithigh, accuracy=1e-8)
    if inithigh - initlow < 0.0
        throw(ArgumentError("inithigh must be larger than initlow!"))
    end

    inits = [x0, xdot0, y0, ydot0]
    t_span = (TMIN, TMAX)
    dlow = initlow
    dhigh = inithigh
    foundlow = false
    foundhigh = false
    ISFROZEN = 10
    frozenloops = 0
    den = 2

    #========== Check if init vals are stable ==============#
    prob = ODEProblem(vectorfield!, inits, t_span, [dlow])
    sol = solve(prob, radau(), reltol=1e-8, verbose=false)

    unstab = false
    for i in 1:length(sol)
        if sol.u[i][3] ≤ 0.0
            unstab = true
        end
    end
    if unstab == false
        foundlow = true
    end

    prob = ODEProblem(vectorfield!, inits, t_span, [dhigh])
    sol = solve(prob, radau(), reltol=1e-8, verbose=false)

    unstab = false
    for i in 1:length(sol)
        if sol.u[i][3] ≤ 0.0
            unstab = true
        end
    end
    if unstab == false
        foundhigh = true
    end

    if foundlow == false && foundhigh == false
        # No stable points in region!
        throw(ArgumentError("No stable diameter in this region!"))
    end

    if foundhigh == true
        return dhigh
    end

    #======================================================#
    println("---------- REDUCTION -------------")
    println("Values converging...")
    println("d_low\t\t\td_high")

    #============= Do reduction ====================#
    # Only happens if low=true, high=false
    while abs(dhigh - dlow) ≥ accuracy
        d = (dhigh + dlow) / den

        if d > dhigh
            # Float error, return low
            return dlow
        end

        prob = ODEProblem(vectorfield!, inits, t_span, [d])
        sol = solve(prob, radau(), reltol=1e-12, verbose=false)
        foundunstab = false
        # Error checking...
        for i in 1:length(sol)
            if sol.u[i][3] ≤ 0.0
                # Unstable d!
                    dhigh = d
                    foundunstab = true
            end
        end

        if foundunstab == false
            dlow = d
        end
        dl = round(dlow, digits=9)
        dh = round(dhigh, digits=9)
        println("$dl\t\t$dh")
    end
    return (dhigh + dlow) / 2
end


function main()
    parsed_args = parse_args(ARGS, s)
    NOPLOT = false
    SCALEDPLOTS = false
    CHOOSE = 0

    if parsed_args["noplot"] != nothing
        NOPLOT = parsed_args["noplot"]
    end
    if parsed_args["scaled-plots"] != nothing
        SCALEDPLOTS = parsed_args["scaled-plots"]
    end
    if parsed_args["choose"] != nothing
        CHOOSE = parsed_args["choose"]
    end

    inits = [x0, xdot0, y0, ydot0]
    t_span = (TMIN, TMAX)


    #====================== Part 1 ====================#
    if CHOOSE == 0 || CHOOSE == 1
        # Container to hold solutions
        solns = []

        println("Running the differential equation with parametric diameter...")

        for d in diameter_vals
            prob = ODEProblem(vectorfield!, inits, t_span, [d])
            sol = solve(prob, radau(), reltol=1e-8, verbose=false)
            push!(solns, sol)
        end

        println("Plotting...")
        plot4(solns, diameter_vals, "time", NOPLOT, SCALEDPLOTS)
        plot4(solns, diameter_vals, "position", NOPLOT, SCALEDPLOTS)

        if CHOOSE == 0
            print("\n\nPress Enter to continue: ")
            readline(stdin)
        end
    end
    #==================================================#


    #====================== Part 2 ====================#
    if CHOOSE == 0 || CHOOSE == 2
        doptim = optimsize(1e-4, 1e-3, 1e-12)
        dop = round(doptim * 1e6, digits=2)
        println("\nCritical diameter: $dop μm")

        proboptim = ODEProblem(vectorfield!, inits, t_span, [doptim])
        soloptim = solve(proboptim, radau(), reltol=1e-8, verbose=false)


        len = length(soloptim)
        tsoptim = soloptim.t
        xsoptim = []
        ysoptim = []
        for i in 1:len
            push!(xsoptim, soloptim.u[i][1])
            push!(ysoptim, soloptim.u[i][3])
        end

        probnoptim = ODEProblem(vectorfield!, inits, t_span, [doptim + 1e-8])
        solnoptim = solve(probnoptim, radau(), reltol=1e-8, verbose=false)

        len = length(solnoptim)
        tsnoptim = solnoptim.t
        xsnoptim = []
        ysnoptim = []
        for i in 1:len
            push!(xsnoptim, solnoptim.u[i][1])
            push!(ysnoptim, solnoptim.u[i][3])
        end

        # Plot critical d and slightly perturbed d
        println("Plotting Critical d...\n")
        plt.clf()
        plt.plot(tsoptim, ysoptim, "r", label="critical d = $dop μm")
        plt.plot(tsnoptim, ysnoptim, "b", label="critical d + 0.01 μm")
        plt.title("Trajectory of Critical Diameter Droplet")
        plt.legend(fontsize=14)
        plt.xlabel("time (s)")
        plt.ylabel("distance (m)")
        plt.ylim(0, 2000)
        plt.show()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
