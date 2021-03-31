module Initialize

include("const.jl")

export init, Model, Grid

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

function init(config) 

    r, hr, ht, u, w, t = 0., 0., 0., 0., 0., 0.
    #  ierr = MPI_Init(argc,argv)

    #Set the cell grid size
    dx = xlen / config["nx_glob"]
    dz = zlen / config["nz_glob"]
    ##############################/
    # BEGIN MPI DUMMY SECTION
    # TODO: (1) GET NUMBER OF MPI RANKS
    #       (2) GET MY MPI RANK ID (RANKS ARE ZERO-BASED INDEX)
    #       (3) COMPUTE MY BEGINNING "I" INDEX (1-based index)
    #       (4) COMPUTE HOW MANY X-DIRECTION CELLS MY RANK HAS
    #       (5) FIND MY LEFT AND RIGHT NEIGHBORING RANK IDs
    ##############################/
    nranks = 1
    myrank = 0
    i_beg = 0
    nx = config["nx_glob"]
    left_rank = 0
    right_rank = 0
    #######################
    # END MPI DUMMY SECTION
    #######################
    println("num_vars = $NUM_VARS")
    #Vertical direction isn't MPI-ized, so the rank's local values = the global values
    k_beg = 0
    nz = config["nz_glob"]
    #masterproc = (myrank == 0)
    data_spec_int = config["data_spec_int"]

    #Allocate the model data
    model = Model(
        zeros(NUM_VARS, nz+2*hs, nx+2*hs),
        zeros(NUM_VARS, nz+2*hs, nx+2*hs),
        zeros(NUM_VARS, nz+1, nx+1),
        zeros(NUM_VARS, nz, nx),
        zeros(nz+2*hs),
        zeros(nz+2*hs),
        zeros(nz+1),
        zeros(nz+1),
        zeros(nz+1),
    )

    grid = Grid(nx, nz, dx, dz)

    #Define the maximum stable time step based on an assumed maximum wind speed
    dt = min(dx,dz) / max_speed * cfl
    #Set initial elapsed model time and output_counter to zero
    etime = 0.
    output_counter = 0.

    #If I'm the master process in MPI, display some grid information
    #  if (masterproc) {
    println("nx_glob, nz_glob: $nx $nz")
    println("dx,dz: $dx $dz")
    println("dt: $dt")
    #  }
    #Want to make sure this info is displayed before further output
    #ierr = MPI_Barrier(MPI_COMM_WORLD)

    # Define quadrature weights and points
    nqpoints = 3
    qpoints = [0.112701665379258311482073460022, 
               0.500000000000000000000000000000, 
               0.887298334620741688517926539980]

    qweights = [0.277777777777777777777777777779,
                0.444444444444444444444444444444,
                0.277777777777777777777777777779]

    #####################################
    # Initialize the cell-averaged fluid state via Gauss-Legendre quadrature
    #####################################
    ########################/
    # TODO: MAKE THESE 2 LOOPS A PARALLEL_FOR
    ########################/
    for k=1:nz+2*hs 
        for i=1:nx+2*hs
            #Use Gauss-Legendre quadrature to initialize a hydrostatic balance + temperature perturbation
            for kk=1:nqpoints
                for ii=1:nqpoints
                    #Compute the x,z location within the global domain based on cell and quadrature index
                    x = (i_beg + i-hs+0.5)*dx + (qpoints[ii]-0.5)*dx
                    z = (k_beg + k-hs+0.5)*dz + (qpoints[kk]-0.5)*dz

                    #Set the fluid state based on the user's specification
                    if data_spec_int == Int(DATA_SPEC_COLLISION)      
                        r, t, u, w, hr, ht = collision(x,z)
                    elseif data_spec_int == Int(DATA_SPEC_THERMAL)
                        r, t, u, w, hr, ht = thermal(x,z)  
                    elseif data_spec_int == Int(DATA_SPEC_MOUNTAIN)
                        r, t, u, w, hr, ht = mountain_waves(x,z) 
                    elseif data_spec_int == Int(DATA_SPEC_TURBULENCE)     
                        r, t, u, w, hr, ht = turbulence(x,z) 
                    elseif data_spec_int == Int(DATA_SPEC_DENSITY_CURRENT) 
                        r, t, u, w, hr, ht = density_current(x,z) 
                    elseif data_spec_int == Int(DATA_SPEC_INJECTION)      
                        r, t, u, w, hr, ht = injection(x,z) 
                    end            
                    #println(r, hr, u, w, ht)
                    #Store into the fluid state array
                    model.state[ID_DENS,k,i] += r                         * qweights[ii]*qweights[kk]
                    model.state[ID_UMOM,k,i] += (r+hr)*u                  * qweights[ii]*qweights[kk]
                    model.state[ID_WMOM,k,i] += (r+hr)*w                  * qweights[ii]*qweights[kk]
                    model.state[ID_RHOT,k,i] += ( (r+hr)*(t+ht) - hr*ht ) * qweights[ii]*qweights[kk]
                end
            end
            for ll=1:NUM_VARS
                model.state_tmp[ll,k,i] = model.state[ll,k,i]
            end
        end
    end
    #Compute the hydrostatic background state over vertical cell averages
    ########################/
    # TODO: MAKE THIS LOOP A PARALLEL_FOR
    ########################/
    for k=1:nz+2*hs
        for kk=1:nqpoints
            z = (k_beg + k-hs+0.5)*dz
            #Set the fluid state based on the user's specification
            if data_spec_int == Int(DATA_SPEC_COLLISION)      
                r, t, u, w, hr, ht = collision(0.,z)
            elseif data_spec_int == Int(DATA_SPEC_THERMAL)
                r, t, u, w, hr, ht = thermal(0.,z)  
            elseif data_spec_int == Int(DATA_SPEC_MOUNTAIN)
                r, t, u, w, hr, ht = mountain_waves(0.,z) 
            elseif data_spec_int == Int(DATA_SPEC_TURBULENCE)     
                r, t, u, w, hr, ht = turbulence(0.,z) 
            elseif data_spec_int == Int(DATA_SPEC_DENSITY_CURRENT) 
                r, t, u, w, hr, ht = density_current(0.,z) 
            elseif data_spec_int == Int(DATA_SPEC_INJECTION)      
                r, t, u, w, hr, ht = injection(0.,z) 
            end
            model.hy_dens_cell[k] = model.hy_dens_cell[k] + hr    * qweights[kk] 
            model.hy_dens_theta_cell[k] = model.hy_dens_theta_cell[k] + hr*ht * qweights[kk]
        end
    end
    #Compute the hydrostatic background state at vertical cell interfaces
    ########################/
    # TODO: MAKE THIS LOOP A PARALLEL_FOR
    ########################/
    for k=1:nz+1
        z = (k_beg + k)*dz
        if data_spec_int == Int(DATA_SPEC_COLLISION)      
            r, t, u, w, hr, ht = collision(0.,z)
        elseif data_spec_int == Int(DATA_SPEC_THERMAL)
            r, t, u, w, hr, ht = thermal(0.,z)  
        elseif data_spec_int == Int(DATA_SPEC_MOUNTAIN)
            r, t, u, w, hr, ht = mountain_waves(0.,z) 
        elseif data_spec_int == Int(DATA_SPEC_TURBULENCE)     
            r, t, u, w, hr, ht = turbulence(0.,z) 
        elseif data_spec_int == Int(DATA_SPEC_DENSITY_CURRENT) 
            r, t, u, w, hr, ht = density_current(0.,z) 
        elseif data_spec_int == Int(DATA_SPEC_INJECTION)      
            r, t, u, w, hr, ht = injection(0.,z) 
        end
        model.hy_dens_int[k] = hr
        model.hy_dens_theta_int[k] = hr * ht
        model.hy_pressure_int[k] = C0 * (hr*ht)^gamma
    end
  
    return model, grid
end 


#This test case is initially balanced but injects fast, cold air from the left boundary near the model top
#x and z are input coordinates at which to sample
#r,u,w,t are output density, u-wind, w-wind, and potential temperature at that location
#hr and ht are output background hydrostatic density and potential temperature at that location
function injection(x , z) 
    hydro_const_theta(z,hr,ht)
    r = 0.
    t = 0.
    u = 0.
    w = 0.
    return r, t, u, w, hr, ht
end


#Initialize a density current (falling cold thermal that propagates along the model bottom)
#x and z are input coordinates at which to sample
#r,u,w,t are output density, u-wind, w-wind, and potential temperature at that location
#hr and ht are output background hydrostatic density and potential temperature at that location
function density_current(x, z) 
    hr, ht = hydro_const_theta(z)
    r = 0.
    t = 0.
    u = 0.
    w = 0.
    t = t + sample_ellipse_cosine(x,z,-20. ,xlen/2,5000.,4000.,2000.)
    return r, t, u, w, hr, ht
end


#x and z are input coordinates at which to sample
#r,u,w,t are output density, u-wind, w-wind, and potential temperature at that location
#hr and ht are output background hydrostatic density and potential temperature at that location
function turbulence(x , z)
    hr, ht = hydro_const_theta(z)
    r = 0.
    t = 0.
    u = 0.
    w = 0.
    # call random_number(u);
    # call random_number(w);
    # u = (u-0.5)*20;
    # w = (w-0.5)*20;
    return r, t, u, w, hr, ht
end


#x and z are input coordinates at which to sample
#r,u,w,t are output density, u-wind, w-wind, and potential temperature at that location
#hr and ht are output background hydrostatic density and potential temperature at that location
function mountain_waves(x, z)
    hr, ht = hydro_const_bvfreq(z, 0.02)
    r = 0.
    t = 0.
    u = 15.
    w = 0.
    return r, t, u, w, hr, ht
end


#Rising thermal
#x and z are input coordinates at which to sample
#r,u,w,t are output density, u-wind, w-wind, and potential temperature at that location
#hr and ht are output background hydrostatic density and potential temperature at that location
function thermal(x, z)
    hr, ht = hydro_const_theta(z)
    r = 0.
    t = 0.
    u = 0.
    w = 0.
    t = t + sample_ellipse_cosine(x,z, 3. ,xlen/2,2000.,2000.,2000.)
    #println(t)
    return r, t, u, w, hr, ht
end


#Colliding thermals
#x and z are input coordinates at which to sample
#r,u,w,t are output density, u-wind, w-wind, and potential temperature at that location
#hr and ht are output background hydrostatic density and potential temperature at that location
function collision!(x, z, r, u, w, t, hr, ht)
    hr, ht = hydro_const_theta!(z)
    r = 0.
    t = 0.
    u = 0.
    w = 0.
    t = t + sample_ellipse_cosine(x,z, 20.,xlen/2,2000.,2000.,2000.)
    t = t + sample_ellipse_cosine(x,z,-20.,xlen/2,8000.,2000.,2000.)
end


#Establish hydrostatic balance using constant potential temperature (thermally neutral atmosphere)
#z is the input coordinate
#r and t are the output background hydrostatic density and potential temperature
function hydro_const_theta(z)
    theta0 = 300.  #Background potential temperature
    exner0 = 1.    #Surface-level Exner pressure
    #Establish hydrostatic balance first using Exner pressure
    ht = theta0                                  #Potential Temperature at z
    exner = exner0 - grav * z / (cp * theta0)   #Exner pressure at z
    p = p0 * exner^(cp/rd)                 #Pressure at z
    rt = (p / C0)^(1. / gamma)              #rho*theta at z
    hr = rt / ht                                  #Density at z
    return hr, ht
end

#Establish hydrstatic balance using constant Brunt-Vaisala frequency
#z is the input coordinate
#bv_freq0 is the constant Brunt-Vaisala frequency
#r and t are the output background hydrostatic density and potential temperature
function hydro_const_bvfreq(z, bv_freq0)
    theta0 = 300.  #Background potential temperature
    exner0 = 1.    #Surface-level Exner pressure
    ht = theta0 * exp(bv_freq0*bv_freq0 / grav * z)                                    #Pot temp at z
    exner = exner0 - grav*grav / (cp * bv_freq0*bv_freq0) * (ht - theta0) / (ht * theta0) #Exner pressure at z
    p = p0 * exner^(cp/rd)                                                         #Pressure at z
    rt = (p / C0)^(1. / gamma)                                                  #rho*theta at z
    hr = rt / ht      
    return hr, ht                                                                    #Density at z
end


#Sample from an ellipse of a specified center, radius, and amplitude at a specified location
#x and z are input coordinates
#amp,x0,z0,xrad,zrad are input amplitude, center, and radius of the ellipse
function sample_ellipse_cosine(x, z, amp, x0, z0, xrad, zrad)
    #Compute distance from bubble center
    dist = sqrt( ((x-x0)/xrad)*((x-x0)/xrad) + ((z-z0)/zrad)*((z-z0)/zrad) ) * pi / 2.
    #If the distance from bubble center is less than the radius, create a cos**2 profile
    if dist <= pi / 2.
        return amp * cos(dist)^2.
    else 
        return 0.
    end
end

end
