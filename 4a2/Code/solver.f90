
      program solver

!     The main body of the CFD solver, calls all subroutines to march towards a
!     flow solution

!     Change to the directory you want to run your case within and execute with
!     "path_to_solver/solver.x < input_casename.txt" to run with screen output
!     or "path_to_solver/solver.x < input_casename.txt > log_casename.txt &" to
!     run in the background

!     Use derived data types defined in a separate module
      use types
      use timestep
      use mesh
      use area
      use check
      use guess
      use loop
      use interp
      use patch
      use secondary
      use writeout

!     Don't use historical implicit variable naming
      implicit none

!     Explicitly declare the required variables
      type(t_appvars) :: av
      type(t_bconds) :: bcs
      type(t_match), allocatable :: p(:)
      type(t_geometry) :: geom
      type(t_grid), allocatable :: g(:) !should it be allocatable
      type(t_grid), allocatable :: g2(:)
      real :: d_max = 1, d_avg = 1
      integer :: nstep, nconv = 5, ncheck = 5, nchange = 5000

!     Read in the data on the run settings
      call read_settings(av,bcs)

!     Determine whether to generate the mesh within this Fortran program or read
!     it directly from a binary file written in Python
      if(av%ni /= -1) then

!         Now the size of the grid is known, the space in memory can be 
!         allocated within the grid type
          call allocate_arrays(av,g,g2,bcs)

!         Read in the case geometry
          call read_geom(av,geom)

!         Set up the mesh coordinates, interpolated between the geometry curves
          call generate_mesh(geom,g,g2)

      else 

!         Read the mesh coordinates directly from file - used for extension
          call read_mesh(av,g,g2,bcs,p)

      end if
      write(*,*) 'done reading mesh'
!     Calculate cell areas and facet lengths
      call calc_areas(g,av)
      call calc_areas(g2,av)

!     Optional output call to inspect the mesh you have generated
      !i guess replace the grid to choose which one to inspect
      call write_output(av,g2,bcs,1)

!     Check that the areas and projected lengths are correct
      call check_mesh(g,av)
      call check_mesh(g2,av)

!     Calculate the initial guess of the flowfield in the domain. There are two
!     options that can be chosen with the input argument "guesstype":
!         1. Uniform flow properties when "guesstype = 1", this is completed
!            for you already, it will allow you to get the solver started but
!            convergence will take more iterations.
!         2. A 1D varying flowfield when "guesstype = 2", assuming isentropic
!            flow in the i-direction allows a calculation of a better
!            approximation to the converged flowfield and so the time to
!            solution will be reduced. You will need to complete this option.
      call flow_guess(av,g,bcs,1)
      

!     Optional output call to inspect the initial guess of the flowfield
      call write_output(av,g,bcs,2)

!     Set the length of the timestep, initially this is a constant based on a 
!     conservative guess of the mach number
      call set_timestep(av,g,bcs)
      
      write(*,*) 'timestep'

!     Open file to store the convergence history. This is human readable during
!     a long run by using "tail -f conv_example.csv" in a terminal window
      open(unit=3,file='conv_' // av%casename // '.csv')

!     Initialise the "stopit" file, during long runs you can request an output
!     is written by setting the value to 1, or to terminate the calculation if
!     set to 2
      open(unit=11,file='stopit')
      write(11,*) 0; close(11);

!     Start the time stepping do loop for "nsteps". This is now the heart of the
!     program, you should aim to program anything inside this loop to operate as
!     efficiently as you can.
      do nstep = 1, av%nsteps

!         Update record of nstep to use in subroutines
          av%nstep = nstep
          if (nstep < 1) then 
                call sub_loop(av, g , bcs, p)
          else if (nstep == 1) then
                call interpolate_mesh(g,g2,bcs,av)
                call set_secondary(av, g2)
                call patch_blocks(g2, p ,av) !in case the interpolation creates discontinuities
                !no longer need to reallocate dt as it is in grid 
                call write_output(av,g2,bcs,3)
            
                call set_timestep(av, g2, bcs)
                call sub_loop(av, g2, bcs, p)

                
          else 
                call sub_loop(av, g2, bcs, p)
          end if 


          !         Stop marching if converged to the desired tolerance "conlim"
          if(d_max < av%d_max .and. d_avg < av%d_avg) then
                write(6,*) 'Calculation converged in', nstep,'iterations'
                exit
          end if


      end do

!     Calculation finished, call "write_output" to write the final, not 
!     necessarily converged flowfield
      write(6,*) 'Calculation completed after', av%nstep,'iterations'
      call write_output(av,g,bcs,3)
!
!     Close open convergence history file
      close(3)

      end program solver


