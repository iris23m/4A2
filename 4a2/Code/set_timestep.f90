      
      subroutine set_timestep(av,g,bcs)

!     This subroutine sets a single value for the time step based on the 
!     stagnation speed of sound and the minimum length scale of any element

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(in) :: g
      type(t_bconds), intent(in) :: bcs
      real :: astag
      real, allocatable :: v_max(:,:)
      real, allocatable :: test(:,:)
      real :: CFL
      integer :: ni, nj

!     Calculate the stagnation speed of sound from the inlet stagnation
!     temperature and gas constants
!     INSERTED
      astag = sqrt(av%gam * av%rgas * bcs%tstag)
      ni = g%ni
      nj = g%nj

!     Assume that the maximum flow speed is also equal to "astag". This will 
!     be pessimistic for subsonic flows but may be optimistic for supersonic 
!     flows. In the latter case the length of the time step as determined by 
!     may need to be reduced by improving this routine or varying the CFL number
!     INSERTED
      allocate(v_max(ni-1,nj-1))
      v_max(1:ni-1, 1:nj-1) = hypot(g%rovx(1:ni-1,1:nj-1)/g%ro(1:ni-1,1:nj-1), g%rovy(1:ni-1,1:nj-1)/g%ro(1:ni-1,1:nj-1))!astag

!     Calculate the timestep using the CFL number and store it in "av%dt"
!     INSERTED

      
      av%dt(1:ni-1,1:nj-1) = av%cfl*g%l_min(1:ni-1,1:nj-1)/(astag+v_max(1:ni-1,1:nj-1))
      !write(*,*) "min sfac" , astag*g%l_min/2


!     Print the calculated timestep and some intermediate values
!     INSERTED
      !write(6,*) 'Calculated timestep: ', av%dt


      end subroutine set_timestep


