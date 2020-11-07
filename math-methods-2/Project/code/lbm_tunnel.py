from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
from math import sqrt
from scipy.interpolate import interp1d
from airfoils import Airfoil

# Define constants:
STEPS_ANIMATION = 10
# Object types: rectangle, circle, horizontal line, vertical line
OBJECT_TYPE = "airfoil"
ANIMATION_TITLE = f"Vorticity Profile of {OBJECT_TYPE}, Lattice Units"
h = 160							# lattice dimensions
w = 350
spacing = 1
height = int(h/spacing)
width = int(w/spacing)
viscosity = 0.018					# fluid viscosity
omega = 1 / (3*viscosity + 0.5)		# "relaxation" parameter
u0 = 0.15						# initial and in-flow speed
four9ths = 4.0/9.0					# abbreviations for lattice-Boltzmann weight factors
one9th   = 1.0/9.0
one36th  = 1.0/36.0

# Initialize all the arrays to steady rightward flow:
n0 = four9ths * (	np.ones((height,width)) - 1.5*u0**2)	# particle densities along 9 directions
nN = one9th * (	np.ones((height,width)) - 1.5*u0**2)
nS = one9th * (	np.ones((height,width)) - 1.5*u0**2)
nE = one9th * (	np.ones((height,width)) + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
nW = one9th * (	np.ones((height,width)) - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
nNE = one36th * (	np.ones((height,width)) + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
nSE = one36th * (	np.ones((height,width)) + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
nNW = one36th * (	np.ones((height,width)) - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
nSW = one36th * (	np.ones((height,width)) - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
rho = n0 + nN + nS + nE + nW + nNE + nSE + nNW + nSW		# macroscopic density
ux = (nE + nNE + nSE - nW - nNW - nSW) / rho				# macroscopic x velocity
uy = (nN + nNE + nNW - nS - nSE - nSW) / rho				# macroscopic y velocity

# Barriers:

def gen_circle(xoff, yoff, rad):
	b = 	np.zeros((height,width), bool)					# True wherever there's a barrier
	for x, y in np.ndindex(b.shape):
		if ((x - xoff)**2 + (y - yoff)**2) < (rad**2):
			b[x, y] = True
		else:
			b[x, y] = False
	return b

def gen_rect(xoff, yoff, xlen, ylen):
	b = 	np.zeros((height,width), bool)					# True wherever there's a barrier
	for x, y in np.ndindex(b.shape):
		if (xoff - xlen/2) <= x <= (xoff + xlen/2) and (yoff - ylen/2) <= y <= (yoff + ylen/2):
			b[x, y] = True
		else:
			b[x, y] = False
	return b

def gen_right_angle(xoff, yoff, xlen, ylen):
	b = np.zeros((height,width), bool)					# True wherever there's a barrier
	for x, y in np.ndindex(b.shape):
		if (xoff - xlen/2) <= x <= (xoff + xlen/2) and (yoff - ylen/2) <= y <= (yoff + ylen/2):
			if x == (xoff - int(xlen/2)) or y == (yoff - int(ylen/2)):
				b[x, y] = True
		else:
			b[x, y] = False
	return b


def gen_airfoil(xoff, yoff, xscale=1, yscale=1):
	b = np.zeros((height,width), bool)

	foil = Airfoil.NACA4('4812', n_points=height)
	xu = np.copy(foil._x_upper)
	yu = np.copy(foil._y_upper)
	xl = np.copy(foil._x_lower)
	yl = np.copy(foil._y_lower)


	xu *= xscale
	xl *= xscale
	yu *= yscale
	yl *= yscale

	# Set origin by shifting
	xu += xoff
	xl += xoff
	yu += yoff
	yl += yoff

	foil._order_data_points()
	plt.plot(xu, yu, xl, yl)


	def y_up(x):
		return interp1d(
			xu, yu, kind="quadratic", bounds_error=False, fill_value="extrapolate"
		)(x)

	def y_low(x):
		return interp1d(
			xl, yl, kind="quadratic", bounds_error=False, fill_value="extrapolate"
		)(x)


	for x,y in np.ndindex(b.shape):
		if float(y_low(y)) <= x <= float(y_up(y)):
			b[x,y] = True
		else:
			b[x,y] = False

	xmin, xmax = b[:,0].min(), b[:,0].max()
	ymin, ymax = b[:,1].min(), b[:,1].max()

	return b, xmin, xmax, ymin, ymax



if OBJECT_TYPE == "circle":
	barrier = gen_circle(int(height/2), int(width/6), 8)
elif OBJECT_TYPE == "rectangle":
	barrier = gen_rect(int(height/2), int(width/3), int(height/10), int(width/8))
elif OBJECT_TYPE == "horizontal line":
	barrier = gen_rect(int(height/2), int(width/3), 1, int(width/8))
elif OBJECT_TYPE == "vertical line":
	barrier = gen_rect(int(height/2), int(width/3), int(height/3), 1)
elif OBJECT_TYPE == "right angle":
	barrier = gen_right_angle(int(height/2), int(width/3), int(height/4), int(width/4))
elif OBJECT_TYPE == "airfoil":
	barrier, xmin, xmax, ymin, ymax = gen_airfoil(20, int(height/2), 180, 140)
	# barrier = gen_airfoil(2, int(height/2), int(width/6), 100, 5, 0.08)




barrierN = 	np.roll(barrier,  1, axis=0)					# sites just north of barriers
barrierS = 	np.roll(barrier, -1, axis=0)					# sites just south of barriers
barrierE = 	np.roll(barrier,  1, axis=1)					# etc.
barrierW = 	np.roll(barrier, -1, axis=1)
barrierNE = 	np.roll(barrierN,  1, axis=1)
barrierNW = 	np.roll(barrierN, -1, axis=1)
barrierSE = 	np.roll(barrierS,  1, axis=1)
barrierSW = 	np.roll(barrierS, -1, axis=1)


# Move all particles by one step along their directions of motion (pbc):
def stream():
	global nN, nS, nE, nW, nNE, nNW, nSE, nSW, barrier
	nN  = 	np.roll(nN,   1, axis=0)	# axis 0 is north-south; + direction is north
	nNE = 	np.roll(nNE,  1, axis=0)
	nNW = 	np.roll(nNW,  1, axis=0)
	nS  = 	np.roll(nS,  -1, axis=0)
	nSE = 	np.roll(nSE, -1, axis=0)
	nSW = 	np.roll(nSW, -1, axis=0)
	nE  = 	np.roll(nE,   1, axis=1)	# axis 1 is east-west; + direction is east
	nNE = 	np.roll(nNE,  1, axis=1)
	nSE = 	np.roll(nSE,  1, axis=1)
	nW  = 	np.roll(nW,  -1, axis=1)
	nNW = 	np.roll(nNW, -1, axis=1)
	nSW = 	np.roll(nSW, -1, axis=1)
	# Use tricky boolean arrays to handle barrier collisions (bounce-back):
	nN[barrierN] = nS[barrier]
	nS[barrierS] = nN[barrier]
	nE[barrierE] = nW[barrier]
	nW[barrierW] = nE[barrier]
	nNE[barrierNE] = nSW[barrier]
	nNW[barrierNW] = nSE[barrier]
	nSE[barrierSE] = nNW[barrier]
	nSW[barrierSW] = nNE[barrier]

# Collide particles within each cell to redistribute velocities (could be optimized a little more):
def collide():
	global rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW
	rho = n0 + nN + nS + nE + nW + nNE + nSE + nNW + nSW
	ux = (nE + nNE + nSE - nW - nNW - nSW) / rho
	uy = (nN + nNE + nNW - nS - nSE - nSW) / rho
	ux2 = ux * ux				# pre-compute terms used repeatedly...
	uy2 = uy * uy
	u2 = ux2 + uy2
	omu215 = 1 - 1.5*u2			# "one minus u2 times 1.5"
	uxuy = ux * uy
	n0 = (1-omega)*n0 + omega * four9ths * rho * omu215
	nN = (1-omega)*nN + omega * one9th * rho * (omu215 + 3*uy + 4.5*uy2)
	nS = (1-omega)*nS + omega * one9th * rho * (omu215 - 3*uy + 4.5*uy2)
	nE = (1-omega)*nE + omega * one9th * rho * (omu215 + 3*ux + 4.5*ux2)
	nW = (1-omega)*nW + omega * one9th * rho * (omu215 - 3*ux + 4.5*ux2)
	nNE = (1-omega)*nNE + omega * one36th * rho * (omu215 + 3*(ux+uy) + 4.5*(u2+2*uxuy))
	nNW = (1-omega)*nNW + omega * one36th * rho * (omu215 + 3*(-ux+uy) + 4.5*(u2-2*uxuy))
	nSE = (1-omega)*nSE + omega * one36th * rho * (omu215 + 3*(ux-uy) + 4.5*(u2-2*uxuy))
	nSW = (1-omega)*nSW + omega * one36th * rho * (omu215 + 3*(-ux-uy) + 4.5*(u2+2*uxuy))
	# Force steady rightward flow at ends (no need to set 0, N, and S components):
	nE[:,0] = one9th * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
	nW[:,0] = one9th * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
	nNE[:,0] = one36th * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
	nSE[:,0] = one36th * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
	nNW[:,0] = one36th * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
	nSW[:,0] = one36th * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)

# Compute curl of the macroscopic velocity field:
def curl(ux, uy):
	return 	np.roll(uy,-1,axis=1) - 	np.roll(uy,1,axis=1) - 	np.roll(ux,-1,axis=0) + 	np.roll(ux,1,axis=0)


def forceAspect(ax,aspect):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

# # Here comes the graphics and animation...
# # theFig = plt.figure(figsize=(5,5))
# theFig = plt.figure()
# # fluidImage = plt.imshow(curl(ux, uy), origin='lower',
# # 	norm=plt.Normalize(-.2,.2), cmap=plt.get_cmap('jet'),
# # 	interpolation='spline16')
# ext = [xmin-10, xmax+10, ymin-0.2, ymax+0.2]
# fluidImage = plt.imshow(curl(ux, uy), origin='lower', cmap=plt.get_cmap('jet'),
# 	interpolation='spline16')#, extent=ext, aspect='auto')
#
#
# bImageArray = np.zeros((height, width, 4), 	np.uint8)	# an RGBA image
# bImageArray[barrier,3] = 255								# set alpha=255 only at barrier sites
# barrierImage = plt.imshow(bImageArray, origin='lower', interpolation='quadric')#, extent=ext, aspect='auto')
# # barrierImage = plt.imshow(bImageArray, origin='lower')#, extent=ext, aspect='auto')
# # forceAspect(ax, 1.0)
#
# plt.title(ANIMATION_TITLE)
# plt.subplots_adjust(top=0.88)
# # plt.xlim(74, 175)
# # plt.ylim(32.85, 33.55)
#
# # Function called for each successive animation frame:
# def nextFrame(arg):
# 	# arg is frame number
#
# 	if arg == 0 or (arg % 15) == 0:
# 		frameName = f"frame{arg}_{OBJECT_TYPE}.png"
# 		plt.savefig(frameName)
#
# 	for step in range(STEPS_ANIMATION):
# 		stream()
# 		collide()
#
# 	fluidImage.set_array(curl(ux, uy))
# 	return (fluidImage, barrierImage)		# return the figure elements to redraw
#
# animate = animation.FuncAnimation(theFig, nextFrame, interval=1, blit=True)
# plt.clim(vmin=0, vmax=1)
# plt.colorbar(fluidImage)
# plt.show()

def julia_plot(f, b, ux, uy, ρ)
