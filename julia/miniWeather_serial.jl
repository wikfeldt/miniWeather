using ArgParse
using NetCDF
include("const.jl")
include("Initialize.jl")
using .Initialize: init, Model, Grid
include("Timestep.jl")
using .Timestep: perform_timestep!

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--nx_glob"
        help = "Total number of cells in the x-direction"
        default = 100
        "--nz_glob"
        help = "Total number of cells in the z-direction"
        default = 50
        "--sim_time"
        help = "How many seconds to run the simulation"
        default = 1000
        "--output_freq"
        help = "How frequently to output data to file (in seconds)"
        default = 10
        "--data_spec_int"
        help = "How to initialize the data"
        default = Int(DATA_SPEC_THERMAL)
    end

    return parse_args(s)
end


#Compute reduced quantities for error checking without resorting to the "ncdiff" tool
function reductions(model, grid)
    mass = 0.0
    te = 0.0
    nx, nz = grid.nx, grid.nz
    dx, dz = grid.dx, grid.dz

    for k = 1:nz
        for i = 1:nx
            r = model.state[i, k, ID_DENS] + model.hy_dens_cell[k]       # Density
            u = model.state[i, k, ID_UMOM] / r                           # U-wind
            w = model.state[i, k, ID_WMOM] / r                           # W-wind
            th = (model.state[i, k, ID_RHOT] + model.hy_dens_theta_cell[k]) / r # Potential Temperature (theta)
            p = C0 * (r * th)^gamma      # Pressure
            t = th / (p0 / p)^(rd / cp)  # Temperature
            ke = r * (u * u + w * w)           # Kinetic Energy
            ie = r * cv * t                # Internal Energy
            mass = mass + r * dx * dz # Accumulate domain mass
            te = te + (ke + r * cv * t) * dx * dz # Accumulate domain total energy
            #println(r, u, w, th)
        end
    end

    #    double glob[2], loc[2];
    #    loc[0] = mass;
    #    loc[1] = te;
    #    int ierr = MPI_Allreduce(loc,glob,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    #    mass = glob[0];
    #    te   = glob[1];
    return mass, te
end

#=
#Output the fluid state (state) to a NetCDF file at a given elapsed model time (etime)
#The file I/O uses parallel-netcdf, the only external library required for this mini-app.
#If it's too cumbersome, you can comment the I/O out, but you'll miss out on some potentially cool graphics
function output(state, etime, grid, nt, ncfile="output.nc")
  integer :: ncid, t_dimid, x_dimid, z_dimid, dens_varid, uwnd_varid, wwnd_varid, theta_varid, t_varid
  integer :: i,k
  integer, save :: num_out = 0
  integer(kind=MPI_OFFSET_KIND) :: len, st1(1),ct1(1),st3(3),ct3(3)
  #If the elapsed time is zero, create the file. Otherwise, open the file
  if etime == 0 
    #Create the file
    nccreate(ncfile, "time", "t", nt, Dict("units"=>"s"))
    nccreate(ncfile, "dens", "x", nx, Dict("units"=>"m"), "z", nz, Dict("units"=>"m"), "time", nt, Dict("units"=>"s"))
    nccreate(ncfile, "uwnd", "x", nx, Dict("units"=>"m"), "z", nz, Dict("units"=>"m"), "time", nt, Dict("units"=>"s"))
    nccreate(ncfile, "wwnd", "x", nx, Dict("units"=>"m"), "z", nz, Dict("units"=>"m"), "time", nt, Dict("units"=>"s"))
    nccreate(ncfile, "theta", "x", nx, Dict("units"=>"m"), "z", nz, Dict("units"=>"m"), "time", nt, Dict("units"=>"s"))    
    dens = zeros(grid.nx, grid.nz, nt)
    uwnd = zeros(grid.nx, grid.nz, nt)
    wwnd = zeros(grid.nx, grid.nz, nt)
    theta = zeros(grid.nx,grid.nz, nt)
  else
    t_old = ncread(ncfile, "t")
    dens_old = ncread(ncfile, "dens")
    uwnd_old = ncread(ncfile, "uwnd")
    wwnd_old = ncread(ncfile, "wwnd")
    theta_old = ncread(ncfile, "theta")
    nt = length(t_old)
    dens = cat(dens_old, dens, dims=3)
    uwnd = cat(uwnd_old, uwnd, dims=3)
    wwnd = cat(wwnd_old, wwnd, dims=3)
    theta = cat(theta_old, theta, dims=3)
    push!(t_old, etime)
  end



  #Store perturbed values in the temp arrays for output
  for k=1:nz
    for i=1:nx
      dens[i,k] = state[i,k,ID_DENS]
      uwnd[i,k] = state[i,k,ID_UMOM] / (hy_dens_cell[k] + state[i,k,ID_DENS])
      wwnd[i,k] = state[i,k,ID_WMOM] / (hy_dens_cell[k] + state[i,k,ID_DENS])
      theta[i,k] = (state[i,k,ID_RHOT] + hy_dens_theta_cell[k]) / (hy_dens_cell[k] + state[i,k,ID_DENS]) - hy_dens_theta_cell[k] / hy_dens_cell[k]
    end
  end

      
end

    #Write the grid data to file with all the processes writing collectively
    ncwrite(dens, ncfile, "dens")
    ncwrite(uwnd, ncfile, "uwnd")
    ncwrite(wwnd, ncfile, "wwnd")
    ncwrite(theta, ncfile, "theta")
    ncwrite(etime, ncfile, "time")
  else


end


  #Write the grid data to file with all the processes writing collectively
  ncwrite(dens, ncfile, "dens")
  ncwrite(uwnd, ncfile, "uwnd")
  ncwrite(wwnd, ncfile, "wwnd")
  ncwrite(theta, ncfile, "theta")
  ncwrite(etime, ncfile, "time")



  #Increment the number of outputs
#  num_out = num_out + 1

  #Deallocate the temp arrays
end 
 =#






function main()

    if !isinteractive()
        config = parse_commandline()
        println(typeof(config))
    else
        config = Dict(
            "nx_glob" => 100,
            "nz_glob" => 50,
            "sim_time" => 1000,
            "output_freq" => 10,
            "data_spec_int" => Int(DATA_SPEC_THERMAL),
        )
    end

    println("Parsed args:")
    for (arg, val) in config
        println("  $arg  =>  $val")
    end

    model, grid = init(config)
    #init(config);
    #println(model)

    mass, te = reductions(model, grid)
    println("mass = $mass, te = $te")

    etime = 0.0
    #output(model, etime)

    time_counter = 0
    direction_switch = true
    for i = 1:grid.nt
        if mod(time_counter, config["output_freq"]) == 0
            println(etime)
        end

        perform_timestep!(model, grid, direction_switch)
        direction_switch = !direction_switch
        etime += grid.dt
        time_counter += 1
    end

end


main()
