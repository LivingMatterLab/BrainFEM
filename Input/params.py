import math
import numpy as np

# Input parameters
nelc = 150;					# Number of elements in circumferential direction
core_frac = .2				# Radial fraction of the core
subcort_frac = .95			# Radial fraction of the cortex
elemDir = 'radial'			# Element directions
nMatProp   = 100;			# Number of materials and properties to create

pert = 0.005            # pertubation of center nodes in cortex (normalized wrt cortical thickness)
pert_width = 0.05       # width of perturbed area (normalized wrt width of plate)

# Properties
E_inner = 1000.;			# Stiffness of inner parts (core and subcortex)
E_ratio = 1.;			# Stiffness ratio between outer (cortex) and inner
kth1_i = 47.0;          	# dtheta1/drho inner (perpendicular to element dir.)
kth1_o = 47.0+0.;          # dtheta1/drho outer (perpendicular to element dir.)
kth2_i = 47.0;          	# dtheta2/drho inner (parallel to element dir.)
kth2_o = 47.0-0.;          # dtheta2/drho outer (parallel to element dir.)
alpha  = 1.656474;				# alpha as in theta = (1+k*rho)^alpha
D = 1100.0;             	# Diffusivity
Grho = 0.121558;          	# Density increase per time
advSpeed = 38.5;  	# Advection speed

# Input for differential Heaviside functions
# H = exp(alpha*(x-center))/(1+ (exp(alpha*(x-center))) )
src_R_center = core_frac					# Density source in r
src_R_alpha  = 5. 

src_T_center = 40;							# Density source in time
src_T_alpha  = 1.5;

adv_R_center = subcort_frac					# Advection in r
adv_R_alpha  = 20. 							

adv_RHO_center = 0.04						# Advection threshold in Rho
adv_RHO_alpha  = 200. 	

# Solver
tEnd = 19.;
tLoad = tEnd;
dt0  = 0.25;
dtMin = 0.00001;
dtMax = dt0;
maxIter = 10;
maxIterInc = 5;
maxStep = 1000;
plotOutput = True
