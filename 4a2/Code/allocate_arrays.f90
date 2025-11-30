      
      subroutine allocate_arrays(av,g,g2,bcs)

!     Allocate memory for all arrays in the grid and bcs datatype, this has been
!     completed for the basic solver. If you add further arrays to the code in
!     the extensions you will need to allocate them here.

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(inout) :: g, g2
      type(t_bconds), intent(inout) :: bcs
      integer :: ni, nj

!     Get the size of the mesh and store locally for convenience
      ni = av%ni; nj = av%nj;
      allocate(g%dt(ni-1,nj-1))

!     Copy the mesh size to the grid datatype
      g%ni = ni; g%nj = nj;
      g2%ni = ni *2-1  
      g2%nj = nj*2-1

!     Wall flag array
      allocate(g%wall(ni,nj))
      allocate(g2%wall(ni*2-1,nj*2-1))
      

!     Arrays to store static conditions at the inlet plane
      allocate(bcs%ro(nj),bcs%p(nj))

!     Node coordinates in the mesh
      allocate(g%x(ni,nj),g%y(ni,nj))
      allocate(g2%x(ni*2-1,nj*2-1),g2%y(ni*2-1,nj*2-1))

!     Cell areas and projected side lengths at the centre of each respectively
      allocate(g%area(ni-1,nj-1),g%lx_i(ni,nj-1),g%ly_i(ni,nj-1), &
          g%lx_j(ni-1,nj),g%ly_j(ni-1,nj), g%l_min(ni-1,nj-1))

      allocate(g2%area(g2%ni-1,g2%nj-1),g2%lx_i(g2%ni,g2%nj-1),g2%ly_i(g2%ni,g2%nj-1), &
          g2%lx_j(g2%ni-1,g2%nj),g2%ly_j(g2%ni-1,g2%nj), g2%l_min(g2%ni-1,g2%nj-1))

!     Primary flow variables in the mesh
      allocate(g%ro(ni,nj),g%rovx(ni,nj),g%rovy(ni,nj),g%roe(ni,nj))

      allocate(g2%ro(ni*2-1,nj*2-1),g2%rovx(ni*2-1,nj*2-1),g2%rovy(ni*2-1,nj*2-1),g2%roe(ni*2-1,nj*2-1))

!     Cell centred primary increments
      allocate(g%dro(ni-1,nj-1),g%drovx(ni-1,nj-1), &
          g%drovy(ni-1,nj-1),g%droe(ni-1,nj-1))

      allocate(g2%dro(g2%ni-1,g2%nj-1),g2%drovx(g2%ni-1,g2%nj-1), &
          g2%drovy(g2%ni-1,g2%nj-1),g2%droe(g2%ni-1,g2%nj-1))

!     Secondary variables stored at the nodes for convenience
      allocate(g%p(ni,nj),g%hstag(ni,nj),g%vx(ni,nj),g%vy(ni,nj))

      allocate(g2%p(ni*2-1,nj*2-1),g2%hstag(ni*2-1,nj*2-1),g2%vx(ni*2-1,nj*2-1),g2%vy(ni*2-1,nj*2-1))


      end subroutine allocate_arrays


