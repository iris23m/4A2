      
      subroutine generate_mesh(geom,g)

!     Create cells of the mesh to cover the domain defined by geometry curves,
!     the values of the node coordinates, x(i,j) and y(i,j) are stored in "g"

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_geometry), intent(in) :: geom
      type(t_grid), intent(inout) :: g
      real :: si_a(geom%ni_a), si_b(geom%ni_b), si(g%ni)
      integer :: ni, nj

!     Declare integers or any extra variables you need here
!     INSERTED
      integer :: i
      integer :: j
      integer :: case, size
      real :: sj(g%nj)
      real :: sum, dsi, dsj, n, centre, ratio, f, df, growthratio, initiallevel, r
      real :: listj(g%nj)
      real, allocatable :: si_1(:)
      real, allocatable :: si_2(:)

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;
      case = 2 !1 is normal mesh density,2 is bend, 3 is bump

!     Calculate the non-dimensional curve lengths of the geometry input and
!     generate linearly spaced points in i-direction at desired mesh resolution
      call dist(geom%x_a,geom%y_a,1,si_a)
      call dist(geom%x_b,geom%y_b,1,si_b)
      if(case == 1 .OR. case==2) then
            call linspace(0.0,1.0,si)
            write(*,*) si
      else if(case == 3) then
            centre = 0.3
            size = NINT(ni * centre)
            write(*,*) size
            allocate(si_1(size))
            allocate(si_2(ni-size))
            n = size -1
            write(*,*) n
            dsi = 1.0/ni * 0.2
            r = 1.2
            
            do 
                  f = dsi*(1 - r**(n))/(1 - r) -centre
                  write(*,*) 'f',  dsi*(1 - r**(n))/(1 - r) -centre
                  df = dsi*( (-n*r**(n-1)*(1 - r) + (1 - r**n)) ) / (1 - r)**2
                  write(*,*) df
                  if (abs(f) < 1e-6) exit
                  r = r - f/df
            end do
            do i = 1,size
                  si_1(i) = centre - dsi * (1-r**(i-1))/(1-r)
            end do
            n = ni - size !-1
            r = 1.2
            do 
                  f = dsi*(1 - r**(n))/(1 - r) -(1-centre)
                  df = dsi*( (-n*r**(n-1)*(1 - r) + (1 - r**n)) ) / (1 - r)**2
                  if (abs(f) < 1e-6) exit
                  r = r - f/df
                  write(*,*) f
            end do
            do i = 1,(ni-size)
                  si_2(i) = dsi * (1-r**(i))/(1-r) + centre
            end do
            si = [si_1(size:1:-1),si_2]
            si(ni) = 1
            si(1) = 0
            write(*,*) si
      end if

!     Interpolate the geometry curves to the required resolution in the 
!     i-direction, this allows the mesh to be refined without altering the 
!     geometry definition file, the data is stored at "j = 1" and "j = nj"
      call interp(si_a,geom%x_a,si,g%x(:,1))
      call interp(si_a,geom%y_a,si,g%y(:,1))
      call interp(si_b,geom%x_b,si,g%x(:,nj))
      call interp(si_b,geom%y_b,si,g%y(:,nj))

!     Calculate the coordinates of all the intermediate points within the mesh.
!     Create a new vector of non-dimensional spacings in the j-direction using 
!     "linspace", loop over the mesh in the i-direction and calculate the
!     intermediate coordinates from a weighted sum of the two boundaries
!     INSERTED
      if(case == 1) then
            call linspace(0.0,1.0,sj)
            write(*,*) sj
      else if(case == 2 .OR. case ==3) then
            dsj = 1.0/nj * 0.2
            r = 1.2
            n = nj
            
            do 
                  f = dsj*(1 - r**(n))/(1 - r) -1.0
                  write(*,*) 'f',  dsj*(1 - r**(n))/(1 - r) -1
                  df = dsj*( (-n*r**(n-1)*(1 - r) + (1 - r**n)) ) / (1 - r)**2
                  write(*,*) df
                  if (abs(f) < 1e-6) exit
                  r = r - f/df
            end do
            do j = 1,nj
                  sj(j) = dsj * (1-r**(j-1))/(1-r)
            end do
            sj(nj) = 1.0
            write(*,*) sj
            
      end if
      do i = 1, ni
            do j = 1, nj 
                  g%x(i,j) = (sj(j)*g%x(i,nj)) + ((1-sj(j))*g%x(i,1))
                  g%y(i,j) = (sj(j)*g%y(i,nj)) + ((1-sj(j))*g%y(i,1))
            end do
      end do  

      


!     In all of the test cases for the basic solver the the "j = 1" and "j = nj"
!     boundaries are walls, for the extensions you may need to return to this
!     and communicate the position of the walls to the solver in a more 
!     flexible way. The "wall" variable is an "ni * nj" logical array where 
!     "true" indicates the node is on a wall.
      g%wall = .false.
      g%wall(:,[1,g%nj]) = .true.

!     Print that the mesh has been created
      write(6,*) 'Interpolated mesh from the bounding geometry curves'
      write(6,*)

      end subroutine generate_mesh


