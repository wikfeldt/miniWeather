using ArgParse
using NetCDF
include("const.jl")
include("Initialize.jl")
using .Initialize: init, Model, Grid
include("Timestep.jl")
using .Timestep: perform_timestep!
using ProgressMeter

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--nx_glob"
        help = "Total number of cells in the x-direction"
        arg_type = Int
        default = 100
        "--nz_glob"
        help = "Total number of cells in the z-direction"
        arg_type = Int
        default = 50
        "--sim_time"
        help = "How many seconds to run the simulation"
        arg_type = Int
        default = 50
        "--output_freq"
        help = "How frequently to output data to file (in seconds)"
        arg_type = Int
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

    for k = 1:nz, i = 1:nx
        r = model.state[hs+i, hs+k, ID_DENS] + model.hy_dens_cell[hs+k]       # Density
        u = model.state[hs+i, hs+k, ID_UMOM] / r                           # U-wind
        w = model.state[hs+i, hs+k, ID_WMOM] / r                           # W-wind
        th = (model.state[hs+i, hs+k, ID_RHOT] + model.hy_dens_theta_cell[hs+k]) / r # Potential Temperature (theta)
        p = C0 * (r * th)^gamma      # Pressure
        t = th / (p0 / p)^(rd / cp)  # Temperature
        ke = r * (u * u + w * w)           # Kinetic Energy
        ie = r * cv * t                # Internal Energy
        mass = mass + r * dx * dz # Accumulate domain mass
        te = te + (ke + ie) * dx * dz # Accumulate domain total energy
        #println(r, u, w, th)
      end

    #    double glob[2], loc[2];
    #    loc[0] = mass;
    #    loc[1] = te;
    #    int ierr = MPI_Allreduce(loc,glob,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    #    mass = glob[0];
    #    te   = glob[1];
    return mass, te
end

function createSnapshot(model, grid)

    snapshot = Array{Float64}(undef, 4, grid.nx, grid.nz)
    #Store perturbed values in the temp arrays for output
    for k = 1:grid.nz
        for i = 1:grid.nx
            snapshot[1, i, k] =  model.state[i, hs+k, ID_DENS]
            snapshot[2, i, k] =
                model.state[i, hs+k, ID_UMOM] /
                (model.hy_dens_cell[hs+k] + model.state[i, hs+k, ID_DENS])
            snapshot[3, i, k] =
                model.state[i, hs+k, ID_WMOM] /
                (model.hy_dens_cell[hs+k] + model.state[i, hs+k, ID_DENS])
            snapshot[4, i, k] =
                (model.state[i, hs+k, ID_RHOT] + model.hy_dens_theta_cell[hs+k]) /
                (model.hy_dens_cell[hs+k] + model.state[i, hs+k, ID_DENS]) -
                model.hy_dens_theta_cell[hs+k] / model.hy_dens_cell[hs+k]
        end
    end

    return snapshot
end


#
#Output the fluid state (state) to a NetCDF file at a given elapsed model time (etime)
#The file I/O uses parallel-netcdf, the only external library required for this mini-app.
#If it's too cumbersome, you can comment the I/O out, but you'll miss out on some potentially cool graphics
function output(snapshots, grid, etimes, ncfile="output.nc")
  #Create new file
  isfile(ncfile) && rm(ncfile)
  _, nx, nz, _ = size(snapshots)
  nt = length(etimes)
  nccreate(ncfile, "dens", "x", collect(range(0,stop=grid.nx*grid.dx, length=grid.nx)), Dict("units"=>"m"), 
                           "z", collect(range(0,stop=grid.nz*grid.dz, length=grid.nz)), Dict("units"=>"m"), 
                           "time", etimes, Dict("units"=>"s"))
  nccreate(ncfile, "uwnd", "x", collect(range(0,stop=grid.nx*grid.dx, length=grid.nx)), Dict("units"=>"m"), 
                           "z", collect(range(0,stop=grid.nz*grid.dz, length=grid.nz)), Dict("units"=>"m"), 
                           "time", etimes, Dict("units"=>"s"))
  nccreate(ncfile, "wwnd", "x", collect(range(0,stop=grid.nx*grid.dx, length=grid.nx)), Dict("units"=>"m"), 
                           "z", collect(range(0,stop=grid.nz*grid.dz, length=grid.nz)), Dict("units"=>"m"), 
                           "time", etimes, Dict("units"=>"s"))
  nccreate(ncfile, "theta", "x", collect(range(0,stop=grid.nx*grid.dx, length=grid.nx)), Dict("units"=>"m"), 
                           "z", collect(range(0,stop=grid.nz*grid.dz, length=grid.nz)), Dict("units"=>"m"), 
                           "time", etimes, Dict("units"=>"s"))

                           #nccreate(ncfile, "uwnd", "x", nx, Dict("units"=>"m"), "z", nz, Dict("units"=>"m"), "time", nt, Dict("units"=>"s"))
  #nccreate(ncfile, "wwnd", "x", nx, Dict("units"=>"m"), "z", nz, Dict("units"=>"m"), "time", nt, Dict("units"=>"s"))
  #nccreate(ncfile, "theta", "x", nx, Dict("units"=>"m"), "z", nz, Dict("units"=>"m"), "time", nt, Dict("units"=>"s"))    
  #exit()
  #Write the grid data to file with all the processes writing collectively
#  ncwrite(etimes, ncfile, "time")
  ncwrite(snapshots[1, :, :, 1:nt], ncfile, "dens")
  ncwrite(snapshots[2, :, :, 1:nt], ncfile, "uwnd")
  ncwrite(snapshots[3, :, :, 1:nt], ncfile, "wwnd")
  ncwrite(snapshots[4, :, :, 1:nt], ncfile, "theta")

end





function main()

    if !isinteractive()
        config = parse_commandline()
    else
        config = Dict(
            "nx_glob" => 100,
            "nz_glob" => 50,
            "sim_time" => 50,
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

    mass0, te0 = reductions(model, grid)
    println("mass = $mass0, te = $te0")

    etime = 0.0
    #output(model, etime)

    #=
    anim_dens = Animation()
    anim_uwnd = Animation()
    anim_wwnd = Animation()
    anim_theta = Animation()
    anim = [anim_dens, anim_uwnd, anim_wwnd, anim_theta]
    #output_gif!(model, etime, grid, anim);
    =#
    output_freq = config["output_freq"]
    sim_time = config["sim_time"]

    nsnaps = floor(Int, sim_time/output_freq) + 1
    snapshots = Array{Float64}(undef, 4, grid.nx, grid.nz, nsnaps)
    etimes = Vector{Float64}(undef, nsnaps)

    # output initial state
    snapshots[:,:,:,1] = createSnapshot(model, grid)


    direction_switch = true
    output_counter = 0.0
    counter = 1
    @showprogress for i in 1:grid.nt
    #while etime < sim_time

        if output_counter >= output_freq
            counter += 1
            output_counter = output_counter - output_freq
            #println(counter, "  ", etime)
            snapshots[:,:,:,counter] = createSnapshot(model, grid)
            etimes[counter] = etime
        end

        perform_timestep!(model, grid, direction_switch)
        direction_switch = !direction_switch
        etime += grid.dt
        output_counter += grid.dt
        #output_gif!(model, etime, grid, anim);

    end

    mass, te = reductions(model, grid)
    println("d_mass: ", (mass - mass0) / mass0)
    println("d_te:   ", (te - te0) / te0)

    output(snapshots, grid, etimes)

end



main()
