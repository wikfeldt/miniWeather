
  #Performs a single dimensionally split time step using a simple low-storate three-stage Runge-Kutta time integrator
  #The dimensional splitting is a second-order-accurate alternating Strang splitting in which the
  #order of directions is alternated each time step.
  #The Runge-Kutta method used here is defined as follows:
  # q*     = q[n] + dt/3 * rhs(q[n])
  # q**    = q[n] + dt/2 * rhs(q*  )
  # q[n+1] = q[n] + dt/1 * rhs(q** )
  function perform_timestep!(model::Model,grid::Grid, direction_switch)
    if (direction_switch) then
      #x-direction first
      model.state_tmp, model.flux, model.tend = semi_discrete_step(model.state, model.state     ,  grid.dt / 3 , DIR_X)
      model.state_tmp, model.flux, model.tend = semi_discrete_step(model.state, model.state_tmp ,  grid.dt / 2 , DIR_X)
      model.state, model.flux, model.tend =  semi_discrete_step(model.state, model.state_tmp ,  grid.dt / 1 , DIR_X)
      #z-direction second
      model.state_tmp, model.flux, model.tend = semi_discrete_step(model.state, model.state    , grid.dt / 3 , DIR_Z)
      model.state_tmp, model.flux, model.tend = semi_discrete_step(model.state, model.state_tmp, grid.dt / 2 , DIR_Z)
      model.state, model.flux, model.tend = semi_discrete_step(model.state, model.state_tmp    , grid.dt / 1 , DIR_Z)
    else
      #z-direction second
      model.state_tmp, model.flux, model.tend = semi_discrete_step(model.state, model.state, grid.dt / 3, DIR_Z)
      model.state_tmp, model.flux, model.tend = semi_discrete_step(model.state, model.state_tmp,  grid.dt / 2, DIR_Z)
      model.state, model.flux, model.tend = semi_discrete_step(model.state, model.state_tmp    , grid.dt / 1, DIR_Z)
      #x-direction first
      model.state_tmp, model.flux, model.tend = semi_discrete_step(model.state, model.state    , grid.dt / 3, DIR_X)
      model.state_tmp, model.flux, model.tend = semi_discrete_step(model.state, model.state_tmp, grid.dt / 2, DIR_X)
      model.state, model.flux, model.tend =     semi_discrete_step(model.state, model.state_tmp, grid.dt / 1, DIR_X)
    end
    return !direction_switch
  end 


  #Perform a single semi-discretized step in time with the form:
  #state_out = state_init + dt * rhs(state_forcing)
  #Meaning the step starts from state_init, computes the rhs using state_forcing, and stores the result in state_out
  function semi_discrete_step(state_init , state_forcing , grid , dir)
    state_out = zeros(1:grid.nx+2*hs,1:grid.nz+2*hs,NUM_VARS)
    flux = zeros(1:grid.nx+1,1:grid.nz+1,NUM_VARS)
    tend = zeros(1:grid.nx,1:grid.nz,NUM_VARS)

    if dir == DIR_X
      #Set the halo values for this MPI task's fluid state in the x-direction
      set_halo_values_x(state_forcing)
      #Compute the time tendencies for the fluid state in the x-direction
      flux, tend = compute_tendencies_x(state_forcing, grid)
    elseif dir == DIR_Z
      #Set the halo values for this MPI task's fluid state in the z-direction
      set_halo_values_z(state_forcing)
      #Compute the time tendencies for the fluid state in the z-direction
      flux, tend = compute_tendencies_z(state_forcing, grid)
    end

    #################################################
    ## TODO: THREAD ME
    #################################################
    #Apply the tendencies to the fluid state
    for ll = 1:NUM_VARS
      for k = 1:grid.nz
        for i = 1:grid.nx
          state_out[i,k,ll] = state_init[i,k,ll] + grid.dt * tend[i,k,ll]
        end
      end
    end
    return state_out, state_forcing, flux, tend
  end 


  #Compute the time tendencies of the fluid state using forcing in the x-direction
  #Since the halos are set in a separate routine, this will not require MPI
  #First, compute the flux vector at each cell interface in the x-direction (including hyperviscosity)
  #Then, compute the tendencies using those fluxes
  function compute_tendencies_x(state, grid)
    flux = zeros(grid.nx+1, grid.nz+1, NUM_VARS)
    tend = zeros(grid.nx, grid.nz, NUM_VARS)
    d3_vals = zeros(NUM_VALS)
    vals = zeros(NUM_VARS)
    stencil = zeros(4)
    #Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * grid.dx / (16*grid.dt)
    #################################################
    ## TODO: THREAD ME
    #################################################
    #Compute fluxes in the x-direction for each cell
    for k = 1:grid.nz
      for i = 1:grid.nx+1
        #Use fourth-order interpolation from four cell averages to compute the value at the interface in question
        for ll = 1:NUM_VARS
          for s = 1:sten_size
            stencil[s] = state[i-hs-1+s,k,ll]
          end
          #Fourth-order-accurate interpolation of the state
          vals[ll] = -stencil[1]/12 + 7*stencil[2]/12 + 7*stencil[3]/12 - stencil[4]/12
          #First-order-accurate interpolation of the third spatial derivative of the state (for artificial viscosity)
          d3_vals[ll] = -stencil[1] + 3*stencil[2] - 3*stencil[3] + stencil[4]
        end

        #Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
        r = vals[ID_DENS] + hy_dens_cell[k]
        u = vals[ID_UMOM] / r
        w = vals[ID_WMOM] / r
        t = ( vals[ID_RHOT] + hy_dens_theta_cell[k] ) / r
        p = C0*(r*t)**gamma

        #Compute the flux vector
        flux[i,k,ID_DENS] = r*u     - hv_coef*d3_vals[ID_DENS]
        flux[i,k,ID_UMOM] = r*u*u+p - hv_coef*d3_vals[ID_UMOM]
        flux[i,k,ID_WMOM] = r*u*w   - hv_coef*d3_vals[ID_WMOM]
        flux[i,k,ID_RHOT] = r*u*t   - hv_coef*d3_vals[ID_RHOT]
      end
    end

    ######
    # TODO: THREAD ME
    #####
    #Use the fluxes to compute tendencies for each cell
    for ll = 1:NUM_VARS
      for k = 1:grid.nz
        for i = 1:grid.nx
          tend[i,k,ll] = -( flux[i+1,k,ll] - flux[i,k,ll] ) / grid.dx
        end
      end
    end
    return flux, tend
  end 


  #Compute the time tendencies of the fluid state using forcing in the z-direction
  #Since the halos are set in a separate routine, this will not require MPI
  #First, compute the flux vector at each cell interface in the z-direction (including hyperviscosity)
  #Then, compute the tendencies using those fluxes
  function compute_tendencies_z(state,grid)
    real(rp), intent(in   ) :: state(1-hs:nx+hs,1-hs:nz+hs,NUM_VARS)
    real(rp), intent(  out) :: flux (nx+1,nz+1,NUM_VARS)
    real(rp), intent(  out) :: tend (nx,nz,NUM_VARS)
    integer :: i,k,ll,s
    real(rp) :: r,u,w,t,p, stencil(4), d3_vals(NUM_VARS), vals(NUM_VARS), hv_coef
    #Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * dz / (16*dt)
    ####
    ## TODO: THREAD ME
    ####
    #Compute fluxes in the x-direction for each cell
    for k = 1:nz+1
      for i = 1:nx
        #Use fourth-order interpolation from four cell averages to compute the value at the interface in question
        for ll = 1:NUM_VARS
          for s = 1:sten_size
            stencil[s] = state[i,k-hs-1+s,ll]
          end
          #Fourth-order-accurate interpolation of the state
          vals[ll] = -stencil[1]/12 + 7*stencil[2]/12 + 7*stencil[3]/12 - stencil[4]/12
          #First-order-accurate interpolation of the third spatial derivative of the state
          d3_vals[ll] = -stencil[1] + 3*stencil[2] - 3*stencil[3] + stencil[4]
        end

        #Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
        r = vals[ID_DENS] + hy_dens_int[k)]
        u = vals[ID_UMOM] / r
        w = vals[ID_WMOM] / r
        t = ( vals(ID_RHOT) + hy_dens_theta_int[k] ) / r
        p = C0*(r*t)**gamma - hy_pressure_int[k]
        #Enforce vertical boundary condition and exact mass conservation
        if k == 1 || k == nz+1
          w                = 0
          d3_vals[ID_DENS] = 0
        end

        #Compute the flux vector with hyperviscosity
        flux[i,k,ID_DENS] = r*w     - hv_coef*d3_vals[ID_DENS]
        flux[i,k,ID_UMOM] = r*w*u   - hv_coef*d3_vals[ID_UMOM]
        flux[i,k,ID_WMOM] = r*w*w+p - hv_coef*d3_vals[ID_WMOM]
        flux[i,k,ID_RHOT] = r*w*t   - hv_coef*d3_vals[ID_RHOT]
      end
    end

    ####
    ## TODO: THREAD ME
    ####
    #Use the fluxes to compute tendencies for each cell
    for ll = 1:NUM_VARS
      for k = 1:nz
        for i = 1:nx
          tend[i,k,ll] = -( flux[i,k+1,ll] - flux[i,k,ll] ) / dz
          if ll == ID_WMOM
            tend[i,k,ID_WMOM] = tend[i,k,ID_WMOM] - state[i,k,ID_DENS]*grav
          end
        end
      end
    end
    return flux, tend
  end 

