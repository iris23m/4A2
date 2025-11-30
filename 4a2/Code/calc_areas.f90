module area
      use types
      implicit none
contains      
      subroutine calc_areas(g, av)

!     Calculate the area of the quadrilateral cells and the lengths of the side
!     facets

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_grid), intent(inout) :: g(:)
      type(t_appvars), intent(inout) :: av
      integer :: ni, nj

!     Declare integers or any extra variables you need here
!     INSERTED
      integer :: i
      integer :: j
      real :: a(2)
      real :: b(2)
      real, dimension(:,:), allocatable :: h_i
      real, dimension(:,:), allocatable :: h_j
      !real, dimension(:,:), allocatable :: min_hi
      !real, dimension(:,:), allocatable :: min_hj
      real :: min_hi
      real :: min_hj
      integer :: n
       
      write(*,*) 'into calc areas'
      do n = 1, av%nn
      !     Get the size of the mesh and store locally for convenience
            ni = g(n)%ni; nj = g(n)%nj;

      !     Calculate the areas of the cells and store in g%area. The area of any
      !     quadrilateral is half of the magnitude of the cross product of the two
      !     vectors that form the diagonals. Check the order of your product so that
      !     the values come out positive! You can do this using two nested loops in
      !     the i and j-directions or in a vectorised way by indexing the coordinate
      !     arrays with lists of indices
      !     INSERTED
            write(*,*) '1'
            do i = 1, ni-1
                  do j = 1, nj-1
                        a(1) = g(n)%x(i+1,j+1) - g(n)%x(i,j)
                        a(2) = g(n)%y(i+1,j+1) - g(n)%y(i,j)
                        b(1) = g(n)%x(i,j+1) - g(n)%x(i+1,j)
                        b(2) = g(n)%y(i,j+1) - g(n)%y(i+1,j)
                        g(n)%area(i,j) = 0.5 * ((a(1)*b(2))-(a(2)*b(1)))
                  end do
            end do
            write(*,*) '1'
      !     Calculate the projected lengths in the x and y-directions on all of the
      !     "i = const" facets and store them in g%lx_i and g%ly_i. When combined
      !     together these two components define a vector that is normal to the facet,
      !     pointing inwards towards the centre of the cell. This is only the case for
      !     the left hand side of the cell, the vector stored in position i,j points
      !     towards the centre of the i,j cell
      !     INSERTED
      !     in this x and y need to swap for it to point normal to the facet
            g(n)%ly_i = -g(n)%x(:,2:nj)+g(n)%x(:, 1:nj-1)
            g(n)%lx_i = g(n)%y(:,2:nj)-g(n)%y(:, 1:nj-1)
            write(*,*) '1'
      !     Now repeat the calculation for the project lengths on the "j=const"
      !     facets. 
      !     INSERTED
            g(n)%ly_j = g(n)%x(2:ni,:)-g(n)%x(1:ni-1, :)
            g(n)%lx_j = -g(n)%y(2:ni,:)+g(n)%y(1:ni-1, :)

      !     Find the minimum length scale in the mesh, this is defined as the length
      !     of the shortest side of all the cells. Call this length "l_min", it is used
      !     to set the timestep from the CFL number. Start by calculating the lengths
      !     of the i and j facets by using the intrinsic function "hypot", this avoids
      !     underflow and overflow errors. Then find the overal minimum value using
      !     both the "min" and "minval" functions.
      !     INSERTED
            write(*,*) '1'
            allocate(h_i(ni,nj-1), h_j(ni-1, nj))
            write(*,*) '1'
            h_i = hypot(g(n)%lx_i, g(n)%ly_i) 
            h_j = hypot(g(n)%lx_j, g(n)%ly_j)  
            do i = 1,ni-1
                  do j = 1,nj-1
                        min_hi = min(h_i(i,j), h_i(i+1,j)) 
                        min_hj = min(h_j(i,j), h_j(i,j+1))
                        g(n)%l_min(i,j) = min(min_hi, min_hj)  !minimum sidelength for each cell
                  end do
            end do

      !     Print the overall minimum length size that has been calculated
            !write(6,*) 'Calculated cell areas and facet lengths'
            !write(6,*) '  Overall minimum element size = ', g(n)%l_min
            !write(6,*)
            deallocate(h_i, h_j)
      end do

      end subroutine calc_areas
end module
