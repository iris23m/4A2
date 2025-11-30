module timestep
      use types
      implicit none
contains
      subroutine set_timestep(av,g,bcs)

!     This subroutine sets a single value for the time step based on the 
!     stagnation speed of sound and the minimum length scale of any element

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(inout) :: g(:)
      type(t_bconds), intent(in) :: bcs
      real :: astag
      real, allocatable :: v_max(:,:)
      real :: CFL
      integer :: ni, nj, n

!     Calculate the stagnation speed of sound from the inlet stagnation
!     temperature and gas constants
!     INSERTED
      astag = sqrt(av%gam * av%rgas * bcs%tstag)
      do n = 1,av%nn
            ni = g(n)%ni
            nj = g(n)%nj

      !     Assume that the maximum flow speed is also equal to "astag". This will 
      !     be pessimistic for subsonic flows but may be optimistic for supersonic 
      !     flows. In the latter case the length of the time step as determined by 
      !     may need to be reduced by improving this routine or varying the CFL number
      !     INSERTED
            allocate(v_max(ni-1,nj-1))
            v_max(1:ni-1, 1:nj-1) = hypot(g(n)%rovx(1:ni-1,1:nj-1)/g(n)%ro(1:ni-1,1:nj-1), g(n)%rovy(1:ni-1,1:nj-1)/g(n)%ro(1:ni-1,1:nj-1))!astag

      !     Calculate the timestep using the CFL number and store it in "av%dt"
      !     INSERTED

            
            g(n)%dt(1:ni-1,1:nj-1) = av%cfl*g(n)%l_min(1:ni-1,1:nj-1)/(astag+v_max(1:ni-1,1:nj-1))
            !write(*,*) "min sfac" , astag*g%l_min/2


      !     Print the calculated timestep and some intermediate values
      !     INSERTED
            !write(6,*) 'Calculated timestep: ', av%dt
            deallocate(v_max)

      end do


      end subroutine set_timestep
end module


