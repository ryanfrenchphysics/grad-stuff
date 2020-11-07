To run:

Assuming you're working on a Linux system, then Julia is available in your official package repositories. Use, i.e., "$ apt install julia" to install.

Before running, you must install the necessary libraries. The easiest way to do this is by opening the REPL by running "$ julia" in your terminal (it's easier to do this when you're in the directory with the Julia file you wish to run). Once in the REPL, you can:

1) press ']' to enter Pkg mode
    - type: "add DifferentialEquations ODEInterfaceDiffEq PyPlot ArgParse"

    - exit Pkg mode by pressing backspace

2) enter: include("HW6.jl")
    - It will take ~1-2 min to precompile (the DifferentialEquations module is enormous). Precompilation only happens once per Julia session, so the next time you run anything with those same modules, there's no overhead.

3) The include should succeed, you'll know if there's no error thrown. You can then enter: main(), which will run the HW6.jl script. This will run it with the default options, but it only takes ~5 seconds to run.

----------------------------------------------------------------

OPTION 2:

1) I would still install the modules as described above (in Pkg mode)

2) From the command line, type:

    julia [ARGS] HW6.jl

    - Args are optional, and they are:

        --noplot, -n:           Stop plots from showing
        --scaled-plots, -s:     Turn on shared-y axes for subplots
        --choose, -c:           Choose to only run one part of the
                                assignment. i.e., -c1 runs only question 1's solution.

3) Again, this will take a while because of precompilation, but I've confirmed that there aren't any issues on Manjaro and OpenSUSE Linux. Let me know if you have any issues!
