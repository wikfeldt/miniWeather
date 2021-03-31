pi        = 3.14159265358979323846264338327   #Pi
grav      = 9.8                               #Gravitational acceleration (m / s^2)
cp        = 1004.0                             #Specific heat of dry air at constant pressure
cv        = 717.0                              #Specific heat of dry air at constant volume
rd        = 287.0                              #Dry air constant for equation of state (P=rho*rd*T)
p0        = 1.e5                              #Standard pressure at the surface in Pascals
C0        = 27.5629410929725921310572974482   #Constant to translate potential temperature into pressure (P=C0*(rho*theta)**gamma)
gamma     = 1.40027894002789400278940027894   #gamma=cp/Rd 
#Define domain and stability-related constants
xlen      = 2.e4    #Length of the domain in the x-direction (meters)
zlen      = 1.e4    #Length of the domain in the z-direction (meters)
hv_beta   = 0.25     #How strong to diffuse the solution: hv_beta \in [0:1]
cfl       = 1.50    #"Courant, Friedrichs, Lewy" number (for numerical stability)
max_speed = 450        #Assumed maximum wave speed during the simulation (speed of sound + speed of wind) (meter / sec)
hs        = 0          #"Halo" size: number of cells beyond the MPI tasks's domain needed for a full "stencil" of information for reconstruction
sten_size = 4          #Size of the stencil used for interpolation

#Parameters for indexing and flags
NUM_VARS = 4           #Number of fluid state variables
ID_DENS  = 1           #index for density ("rho")
ID_UMOM  = 2           #index for momentum in the x-direction ("rho * u")
ID_WMOM  = 3           #index for momentum in the z-direction ("rho * w")
ID_RHOT  = 4           #index for density * potential temperature ("rho * theta")
DIR_X = 1              #Integer constant to express that this operation is in the x-direction
DIR_Z = 2              #Integer constant to express that this operation is in the z-direction

@enum DATA_SPEC begin
    DATA_SPEC_COLLISION       = 1
    DATA_SPEC_THERMAL         = 2
    DATA_SPEC_MOUNTAIN        = 3
    DATA_SPEC_TURBULENCE      = 4
    DATA_SPEC_DENSITY_CURRENT = 5
    DATA_SPEC_INJECTION       = 6
end

mutable struct Model
    state::Array{Float64,3}             
    state_tmp::Array{Float64,3}         
    flux::Array{Float64,3}              
    tend::Array{Float64,3}              
    hy_dens_cell::Vector{Float64}       
    hy_dens_theta_cell::Vector{Float64} 
    hy_dens_int::Vector{Float64}        
    hy_dens_theta_int::Vector{Float64}  
    hy_pressure_int::Vector{Float64}
end

struct Grid
    nx::Int64
    nz::Int64
    dx::Int64
    dz::Int64
end