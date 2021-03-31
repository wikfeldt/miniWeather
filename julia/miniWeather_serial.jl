using ArgParse
include("const.jl")
include("Initialize.jl")

using .Initialize: init, Model, Grid

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
function reductions(model::Model, grid::Grid) 
  mass = 0.0
  te   = 0.0
  nx, nz = grid.nx, grid.nz
  dx, dz = grid.dx, grid.dz 

  for k = 1:nz
      for i = 1:nx
            r  =   model.state[ID_DENS,k,i] + model.hy_dens_cell[k]       # Density
            u  =   model.state[ID_UMOM,k,i] / r                           # U-wind
            w  =   model.state[ID_WMOM,k,i] / r                           # W-wind
            th = (model.state[ID_RHOT, k,i] + model.hy_dens_theta_cell[k]) / r # Potential Temperature (theta)
            p  = C0*(r*th)^gamma      # Pressure
            t  = th / (p0/p)^(rd/cp)  # Temperature
            ke = r*(u*u+w*w)           # Kinetic Energy
            ie = r*cv*t                # Internal Energy
            mass = mass + r            *dx*dz # Accumulate domain mass
            te   = te   + (ke + r*cv*t)*dx*dz # Accumulate domain total energy
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
function output(state, etime, grid)
  using NetCDF
  integer :: ncid, t_dimid, x_dimid, z_dimid, dens_varid, uwnd_varid, wwnd_varid, theta_varid, t_varid
  integer :: i,k
  integer, save :: num_out = 0
  integer(kind=MPI_OFFSET_KIND) :: len, st1(1),ct1(1),st3(3),ct3(3)
  !Temporary arrays to hold density, u-wind, w-wind, and potential temperature (theta)
  real(rp), allocatable :: dens(:,:), uwnd(:,:), wwnd(:,:), theta(:,:)
  real(rp) :: etimearr(1)
  #Allocate some (big) temp arrays
  dens = zeros(grid.nx, grid.nz)
  uwnd = xeros(grid.nx, grid.nz)
  wwnd = zeros(grid.nx, grid.nz)
  theta = zeros(grid.nx,grid.nz)

  #If the elapsed time is zero, create the file. Otherwise, open the file
  if etime == 0 
    !Create the file
    nccreate("output.nc", varname, "t", collect(11:20), "t", 20, Dict("units"=>"s"), atts=attribs)

    call ncwrap( nf90mpi_create( MPI_COMM_WORLD , 'output.nc' , nf90_clobber , MPI_INFO_NULL , ncid ) , __LINE__ )
    !Create the dimensions
    len=nf90_unlimited; call ncwrap( nfmpi_def_dim( ncid , 't' , len , t_dimid ) , __LINE__ )
    len=nx_glob       ; call ncwrap( nfmpi_def_dim( ncid , 'x' , len , x_dimid ) , __LINE__ )
    len=nz_glob       ; call ncwrap( nfmpi_def_dim( ncid , 'z' , len , z_dimid ) , __LINE__ )
    !Create the variables
    call ncwrap( nfmpi_def_var( ncid , 't' , nf90_double , 1 , (/ t_dimid /) , t_varid ) , __LINE__ )
    call ncwrap( nfmpi_def_var( ncid , 'dens'  , nf90_double , 3 , (/ x_dimid , z_dimid , t_dimid /) ,  dens_varid ) , __LINE__ )
    call ncwrap( nfmpi_def_var( ncid , 'uwnd'  , nf90_double , 3 , (/ x_dimid , z_dimid , t_dimid /) ,  uwnd_varid ) , __LINE__ )
    call ncwrap( nfmpi_def_var( ncid , 'wwnd'  , nf90_double , 3 , (/ x_dimid , z_dimid , t_dimid /) ,  wwnd_varid ) , __LINE__ )
    call ncwrap( nfmpi_def_var( ncid , 'theta' , nf90_double , 3 , (/ x_dimid , z_dimid , t_dimid /) , theta_varid ) , __LINE__ )
    !End "define" mode
    call ncwrap( nfmpi_enddef( ncid ) , __LINE__ )
  else
    !Open the file
    call ncwrap( nfmpi_open( MPI_COMM_WORLD , 'output.nc' , nf90_write , MPI_INFO_NULL , ncid ) , __LINE__ )
    !Get the variable IDs
    call ncwrap( nfmpi_inq_varid( ncid , 'dens'  ,  dens_varid ) , __LINE__ )
    call ncwrap( nfmpi_inq_varid( ncid , 'uwnd'  ,  uwnd_varid ) , __LINE__ )
    call ncwrap( nfmpi_inq_varid( ncid , 'wwnd'  ,  wwnd_varid ) , __LINE__ )
    call ncwrap( nfmpi_inq_varid( ncid , 'theta' , theta_varid ) , __LINE__ )
    call ncwrap( nfmpi_inq_varid( ncid , 't'     ,     t_varid ) , __LINE__ )
  endif

  !Store perturbed values in the temp arrays for output
  do k = 1 , nz
    do i = 1 , nx
      dens (i,k) = state(i,k,ID_DENS)
      uwnd (i,k) = state(i,k,ID_UMOM) / ( hy_dens_cell(k) + state(i,k,ID_DENS) )
      wwnd (i,k) = state(i,k,ID_WMOM) / ( hy_dens_cell(k) + state(i,k,ID_DENS) )
      theta(i,k) = ( state(i,k,ID_RHOT) + hy_dens_theta_cell(k) ) / ( hy_dens_cell(k) + state(i,k,ID_DENS) ) - hy_dens_theta_cell(k) / hy_dens_cell(k)
    enddo
  enddo

  !Write the grid data to file with all the processes writing collectively
  st3=(/i_beg,k_beg,num_out+1/); ct3=(/nx,nz,1/); call ncwrap( nfmpi_put_vara_double_all( ncid ,  dens_varid , st3 , ct3 , dens  ) , __LINE__ )
  st3=(/i_beg,k_beg,num_out+1/); ct3=(/nx,nz,1/); call ncwrap( nfmpi_put_vara_double_all( ncid ,  uwnd_varid , st3 , ct3 , uwnd  ) , __LINE__ )
  st3=(/i_beg,k_beg,num_out+1/); ct3=(/nx,nz,1/); call ncwrap( nfmpi_put_vara_double_all( ncid ,  wwnd_varid , st3 , ct3 , wwnd  ) , __LINE__ )
  st3=(/i_beg,k_beg,num_out+1/); ct3=(/nx,nz,1/); call ncwrap( nfmpi_put_vara_double_all( ncid , theta_varid , st3 , ct3 , theta ) , __LINE__ )

  !Only the master process needs to write the elapsed time
  !Begin "independent" write mode
  call ncwrap( nfmpi_begin_indep_data(ncid) , __LINE__ )
  !write elapsed time to file
  if (masterproc) then
    st1=(/num_out+1/); ct1=(/1/); etimearr(1) = etime; call ncwrap( nfmpi_put_vara_double( ncid , t_varid , st1 , ct1 , etimearr ) , __LINE__ )
  endif
  !End "independent" write mode
  call ncwrap( nfmpi_end_indep_data(ncid) , __LINE__ )

  !Close the file
  call ncwrap( nf90mpi_close(ncid) , __LINE__ )

  !Increment the number of outputs
  num_out = num_out + 1

  !Deallocate the temp arrays
  deallocate(dens )
  deallocate(uwnd )
  deallocate(wwnd )
  deallocate(theta)
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
            "output" => 10,
            "data_spec_int" => Int(DATA_SPEC_THERMAL)
            )
    end

    println("Parsed args:")
    for (arg,val) in config
        println("  $arg  =>  $val")
    end

    model, grid = init(config)
    #init(config);
    #println(model)

    mass, te = reductions(model, grid)
    println("mass = $mass, te = $te")
end


main()