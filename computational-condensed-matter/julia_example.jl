# Comments in Julia begin with a pound sign
#=
    Multiline comments
=#

# Use the module SpecialFunctions
using SpecialFunctions

#=
    One awesome thing about Julia is that it supports unicode characters using LaTeX notation.
    For example, if one wants to use the constant pi, we can write '\pi' and then hit tab
    to insert. We can also use these sorts of things when naming variables/functions:
    Example -   radius = 5
                radius² = radius * radius
                --> This is done by 'radius\^2' + TAB
=#


#=
    Julia is much more like C than Python or C++ as it does not allow for OOP.
    We can, however, create structures that hold the necessary information that we need.
=#

# Create struct to hold integral limits and number of steps
mutable struct IntegralInfo
    # While in Python we had no control over datatypes of these variables (we had to use a try/except block),
    # We can ensure type safety (which, in turn, increases speed). Low and high must be real numbers, and Julia
    # has a type called Real just for that. Real includes both int and float vals.
    low::Real
    high::Real
    steps::Int64            # Steps is 64-bit integer
    # However, we also need steps to be a positive integer. So we will use an inline function to enforce this:

    #checkifeven(steps) = if (steps % 2 != 0) steps += 1 end

    ### IDK Why this inline is throwing a method error?
end

# ^^^ All functions, structs, loops, etc. must be finished with an 'end' statement

# Create our Simpsons Rule Integration function. Our first argument is a function and our second is the integral's info.
function simpsonsintegrate(func::Function, info::IntegralInfo)
    Δx = (info.high - info.low) / info.steps

    # Create empty arrays for x and f(x), of type Float64
    x = Float64[]
    f = Float64[]

    for i in 1:info.steps + 1
        # To modify values of an object in-place, we place a ! after the function name
        # So, instead of x[i] = value, we can do push!(x, value)
        push!(x, info.low + ((i - 1) * Δx))
        push!(f, func(x[i]))
    end

    result::Float64 = 0.0
    for i in 1:info.steps + 1
        # '||' is boolean or, and '&&' is boolean and, like in C and C++.
        if (i == 1 || i == info.steps + 1)
            result += f[i]
        elseif (i % 2 == 0)
            result += 4 * f[i]
        else
            result += 2 * f[i]
        end
    end
    result = result * (Δx / 3)

    # In Julia, the last value calculated is the returned value. So we don't need 'return result' here,
    # because result would be returned automatically. But this is up to the programmer. I tend to use return
    # because it is easier for other programmers to understand what's going on.
    return result
end

#=
    Let's create some functions to integrate. Another awesome thing about Julia is that most special functions
    commonly used in physics, like cosh, are already included. If they aren't, they are included in the module
    SpecialFunctions
=#
function func1(x::Float64)
    # Unlike most languages, Julia uses ^ for exponentiation
    x^2
end

function func2(x::Float64)
    cosh(x)
end


# Create main function
function main()
    info1 = IntegralInfo(0, 5, 100)
    info2 = IntegralInfo(0, 10, 500)

    result1 = simpsonsintegrate(func1, info1)
    result2 = simpsonsintegrate(func2, info2)

    # Formatted printing in Julia uses $. However, one can also use the macro @printf to use the same
    # type of printf that C uses.

    # Also, it's important to know that Julia uses single quotes for characters and double quotes for strings

    println("Integral 1: x^2 from 0 to 5: $result1")
    println("Integral 2: cosh(x) from 0 to 10: $result2")

    # If we were to pass anything other than real numbers to simpsonsintegrate, the program would crash.

    # Empty return statement so that we aren't returning any nonsense
    return
end


# Finally, let's call our main function
main()



#=
    The workflow for Julia is quite different from other languages. To run this program, open Julia's REPL
    by typing 'julia' into the command line. This brings up the REPL. Now, to go to shell mode (so we can
    navigate to the location of our file) by pressing semicolon. Navigate to our file with cd.
    Finally, we run our program by typing include("julia_example.jl").
=#
