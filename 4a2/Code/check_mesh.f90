 module check
      use types
      implicit none
contains     
      subroutine check_mesh(g,av)

!     Check the cell area and facet length calculations before attempting to
!     solve the flow, make sure you do this for both the "bump" and "bend" test
!     cases

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_grid), intent(inout) :: g(:)
      type(t_appvars), intent(inout) :: av
      real :: area_min, dx_error, dy_error, tol
      integer :: ni, nj, n
      real :: tmp

      do n = 1, av%nn
      !     Get the size of the mesh and store locally for convenience
            ni = g(n)%ni; nj = g(n)%nj;

      !     Exact checking of floating point numbers never goes well, define a
      !     small tolerance to use for a comparative operation instead
            tmp = minval(g(n)%l_min)
            tol = 1e-4 * tmp

      !     Check that all of the cell areas are positive, either with the intrinsic
      !     "minval" function or with nested do loops. Print the output to the screen
      !     and flag negative numbers as an error with an if statement to "stop" the
      !     program
      !     INSERTED
            area_min = minval(g(n)%area)
            write(6,*) 'minimum cell area: ', area_min
            if (area_min <= 0) then
                  stop
            end if

      !     Next check that the sum of the edge vectors around every quadrilateral is 
      !     very nearly zero in both the x and y-coordinate directions. You can
      !     complete this with some elementwise addition of the arrays and use of the
      !     "maxval" and "abs" intrinsic functions.
      !     INSERTED
            dx_error = maxval(-g(n)%lx_i(1:ni-1,1:nj-1)- g(n)%lx_j(1:ni-1,1:nj-1) +g(n)%lx_i(2:ni,1:nj-1) + g(n)%lx_j(1:ni-1,2:nj))
            dy_error = maxval(-g(n)%ly_i(1:ni-1,1:nj-1)- g(n)%ly_j(1:ni-1,1:nj-1) +g(n)%ly_i(2:ni,1:nj-1) + g(n)%ly_j(1:ni-1,2:nj))
            
            if (abs(dx_error)>tol) then
                  write(6,*) 'facet vectors dont sum to 0 in x' , maxval(g(n)%lx_i(1:ni-1,1:nj-1)+ g(n)%lx_j(1:ni-1,1:nj-1) +g(n)%lx_i(2:ni,1:nj-1) + g(n)%lx_j(1:ni-1,2:nj))
            else if (abs(dy_error)>tol) then
                  write(6,*) 'facet vectors dont sum to 0 in y'
            else 
                  write(6,*) 'Check complete: facet vectors sum to 0'
            end if


      !     It may be worthwhile to complete some other checks, the prevous call to
      !     the "write_output" subroutine has written a file that you can read and
      !     postprocess using the Python script plot_mesh.py. This program also has
      !     access to all of the mesh parameters used within the solver that you could
      !     inspect graphically.

      !     Print a blank line
            write(6,*)
      end do

      end subroutine check_mesh
end module
