module loop
      use types
      implicit none
contains
subroutine sub_loop(av,g,bcs,p)
      use types 
      use secondary
      use bconds
      use euler
      use patch
      use timestep
      use conv
      use checkstop
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(inout) :: g(:)
      type(t_bconds), intent(inout) :: bcs
      type(t_match), intent(inout) :: p(:)
      real :: d_max = 1, d_avg = 1
      integer :: nstep, nconv = 1, ncheck = 5, n
      
      if (mod(nstep,10)==0) then
            call set_timestep(av,g,bcs)
      end if
      
      
!         Calculate secondary flow variables used in conservation equations
      call set_secondary(av,g)

!         Apply inlet and outlet values at the boundaries of the domain
      call apply_bconds(av,g,bcs)
      
!         Perform the timestep to update the primary flow variables
      call euler_iteration(av,g)

      call patch_blocks(g, p, av)
      
!         Write out summary every "nconv" steps and update "davg" and "dmax" 
      if(mod(av%nstep,nconv) == 0) then
            
            call check_conv(av,g,d_avg,d_max)
            
      end if

!         Check the solution hasn't diverged or a stop has been requested 
!         every "ncheck" steps
      if(mod(av%nstep,ncheck) == 0) then
            call check_stop(av,g)
      end if


end subroutine sub_loop

end module
