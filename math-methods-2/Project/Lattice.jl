module Lattice
#=
References shorthand:
    Sukrop, Thorne - Lattice Boltzmann Modeling: ST
    Kruger, Kusumaatmaja - Lattice Boltzmann Method: KK
=#

#=
    D2Q9 Geometry Directions (ST p32):

        6   2   5

        3   0   1

        7   4   8

    The velocity magnitude of directions 1-4 is 1 lattice step per time step ,1 lu ts⁻¹, and in directions 5-8 is √2 lu ts⁻¹

    x,y convention of velocities:

        (-1,1)  (0,1)   (1,1)

        (-1,0)  (0,0)   (1,0)

        (-1,-1) (0,-1)  (1,-1)
=#

# Directions: Center, 1-4, 5-8
directions = [
    [0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
    [1, 1], [-1, 1], [-1, -1], [1, -1]
]
export directions

#=
(ST p35)
Collision of the fluid particles is considered as a relaxation towards a local equilibrium and the D2Q9 equilibrium distribution function f eq is de-
fined as:

    fₐ(x) = wₐ ρ(x){
    3(eₐ*u)/c² + (9/2)(eₐ*u)²/c⁴ - (3/2)(u²/c²)
    }

where the weights wa are 4/9 for the rest particles (a = 0), 1/9 for a = 1, 2,
3, 4, and 1/36 for a = 5, 6, 7, 8, and c is the basic speed on the lattice (1 lu
ts-1 in the simplest implementation). Note that if the macroscopic velocity
u = 0, the equilibrium fa are simply the weights times the fluid density.

Weights:

    (1/36)  (1/9)   (1/36)

    (1/9)   (4/9)   (1/9)

    (1/36)  (1/9)   (1/36)
=#

weights = [
    4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
]
export weights


#=
(ST p44)
Boundary conditions:
    Bounceback:
        6 2 5       8 4 7
        3 0 1   ->  1 0 3
        7 4 8       5 2 6

    Bounceback y, slip x:
        6 2 5       7 4 8
        3 0 1   ->  3 0 1
        7 4 8       6 2 5

    Bouncback x, slip y:
        6 2 5       5 2 6
        3 0 1   ->  1 0 3
        7 4 8       8 4 7

Note all arrays go center, 1-4, 5-8
=#
bc_bounceback = [0, 3, 4, 1, 2, 7, 8, 5, 6]
bc_bounceyslipx = [0, 1, 4, 3, 2, 8, 7, 6, 5]
bc_bouncexslipy = [0, 3, 2, 1, 4, 6, 5, 8, 7]
export bc_bounceback, bc_bouceyslipx, bc_bouncexslipy

#=
For future BCs, create arrays holding values for individual columns and rows, i.e.:

                row1 = [2, 5, 6]
                row2 = [0, 1, 3]
    6 2 5       row3 = [7, 8, 4]
    3 0 1   ->  col1 = [3, 6, 7]
    7 4 8       col2 = [0, 2, 4]
                col3 = [1, 5, 8]
=#

row1 = [2, 5, 6]
row2 = [0, 1, 3]
row3 = [7, 8, 4]
col1 = [3, 6, 7]
col2 = [0, 2, 4]
col3 = [1, 5, 8]
export row1, row2, row3, col1, col2, col3

end # module Lattice
