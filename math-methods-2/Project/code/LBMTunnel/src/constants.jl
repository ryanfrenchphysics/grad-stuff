#========= CONSTANTS =============#
# Channel width and height
CHEIGHT = 100 # m
CWIDTH = 200 # m
DIAMETER = CHEIGHT # m/s

# Lattice height, width
δx = 350
δy = convert(Int64, δx * 0.5)

# Speed of sound
c = 1.0 / 3.0

# Initial in-flow speed (Mach #)
u₀ = 0.3

# Kinematic viscosity
ν = 0.02

# Reynold's number (initial)
Re = (u₀ * DIAMETER) / (ν)

# Relaxation parameter
ω = 1.0 / (3.0 * ν + 0.5)

# Length conversion
LAT_X_TO_SI = δx / CWIDTH
LAT_Y_TO_SI = δy / CHEIGHT


# Weight factors
four9ths = 4.0 / 9.0
one9th = 1.0 / 9.0
one36th = 1.0 / 36.0

# BC types
BC_BOUNCEBACK = 1


STEPS_ANIMATION = 15
