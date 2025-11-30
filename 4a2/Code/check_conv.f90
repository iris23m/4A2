module conv
      use types
      use routines
      implicit none 

contains
      subroutine check_conv(av,g,d_avg,d_max)

!     This subroutine checks the residuals in all primary flow variables and
!     prints their values, you should not need to change this subroutine

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(in) :: g(:)
      real, intent(out) :: d_avg, d_max
      real, allocatable :: dro(:,:), droe(:,:), drovx(:,:), drovy(:,:)
      integer :: ijx_max(2), ijy_max(2), ij_max(2), ncells, n
      real :: dro_max, drovx_max, drovy_max, droe_max, dro_avg, drovx_avg, &
          drovy_avg, droe_avg, flow_ratio
      character(len=100) :: fmt_step

      do n = 1,av%nn

            allocate(dro(g(n)%ni-1,g(n)%nj-1), droe(g(n)%ni-1,g(n)%nj-1), drovx(g(n)%ni-1,g(n)%nj-1), drovy(g(n)%ni-1,g(n)%nj-1))
      !     Get the number of cells from the size of the residual arrays
            ncells = size(g(n)%dro)

      !     Use "abs" to make all residual values positive and store locally
            dro = abs(g(n)%dro); droe = abs(g(n)%droe);
            drovx = abs(g(n)%drovx); drovy = abs(g(n)%drovy);

      !     Calculate the mean changes for each variable
            dro_avg = sum(abs(dro)) / (ncells * av%ro_ref)
            droe_avg = sum(abs(droe)) / (ncells * av%roe_ref)
            drovx_avg = sum(abs(drovx)) / (ncells * av%rov_ref)
            drovy_avg = sum(abs(drovy)) / (ncells * av%rov_ref)

      !     Find the maximum value of change for the momenta and the positions
            dro_max = maxval(dro) / av%ro_ref; droe_max = maxval(droe) / av%roe_ref;
            ijx_max = maxloc(drovx); ijy_max = maxloc(drovy);
            drovx_max = drovx(ijx_max(1),ijx_max(2)) / av%rov_ref
            drovy_max = drovy(ijy_max(1),ijy_max(2)) / av%rov_ref

      !     Store single values as the maximum of either the x or y-momentum
            if(drovx_avg > drovy_avg) then
            d_max = drovx_max; d_avg = drovx_avg; ij_max = ijx_max;
            else
            d_max = drovy_max; d_avg = drovy_avg; ij_max = ijy_max;
            end if

      !     Write the average and maximum changes in the primary variables to unit 3
      !     for convergenge plotting
            write(3,'(i13,8e15.6)') av%nstep, dro_avg, droe_avg, drovx_avg, &
            drovy_avg, dro_max, droe_max, drovx_max, drovy_max

      !     Write a short human readable output summary to the screen.
            write(6,*) 'Time step number ', av%nstep, ' Block, ' , n
            fmt_step = '(a,e10.3,a,i4,a,i4,a,e10.3)'
            write(*,fmt_step) '   d_max =', d_max, ' at i =', ij_max(1), ', j =', &
            ij_max(2), ', d_avg =', d_avg

            deallocate(dro, droe, drovx, drovy)
      end do

      end subroutine check_conv

end module


