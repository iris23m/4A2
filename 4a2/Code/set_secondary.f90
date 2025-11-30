      
module secondary
      use types
      implicit none
contains
      subroutine set_secondary(av,g)

!     This subroutine calculates the secondary flow variables from the primary
!     ones at every node

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g(:)
      real, dimension(:,:), allocatable :: T, M_squared
      integer n 
      write(*,*) 'setsec', allocated(g(1)%rovx)
      do n =1,av%nn
      !     Define any further variables you may need
      !     INSERTED
            

      !     The primary flow variables are "ro", "roe", "rovx" and "rovy", these are 
      !     the conserved quantities within the Euler equations. Write the code to
      !     calculate the secondary flow variables which are the velocity components
      !     "vx" and "vy", the static pressure "p" and the stagnation enthalpy
      !     "hstag". These are needed at every timestep, there is no need for any 
      !     loops as the operations can be performed elementwise, although you may
      !     wish to define some intermediate variables to improve readability.
      !     INSERTED
            

            g(n)%vx = g(n)%rovx/g(n)%ro
            g(n)%vy = g(n)%rovy/g(n)%ro
            T = (g(n)%roe - 0.5*g(n)%ro*(g(n)%vx**2 + g(n)%vy**2))/(av%cv*g(n)%ro)
            g(n)%p = g(n)%ro * av%rgas * T
            M_squared = (g(n)%vx**2+ g(n)%vy**2)/(av%gam*av%rgas*T)
            g(n)%hstag = av%cp *(1+0.5*(av%gam-1)*M_squared)*T

            deallocate(T,M_squared)
      end do

      end subroutine set_secondary
end module


