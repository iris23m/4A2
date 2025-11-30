 
 module mesh
      use types
      implicit none
contains
      subroutine read_mesh(av,g,g2,bcs,p)

!     Read a multi-block mesh directly from a binary file, you should not need
!     to change this subroutine

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), allocatable, intent(out) :: g(:),g2(:)
      type(t_bconds), intent(inout) :: bcs
      type(t_match), allocatable, intent(out) :: p(:)
      integer :: n, ni, nj, m, i ,j 

!     Open the file and assign to unit 2
      open(2,file='mesh_' // av%casename // '.bin',form='unformatted', &
          access='stream')

!     Read the number of blocks and allocate the size of "g" 
      read(2) av%nn; allocate(g(av%nn)); allocate(g2(av%nn));
      write(6,*) 'Read multi-block mesh from file'
      write(6,*) '  Number of blocks nn =', av%nn

      

!     Loop over every block and read the mesh coordinates
      do n = 1,av%nn

            

!         Read the size of the block 
          read(2) ni, nj
          g(n)%ni = ni; g(n)%nj = nj;
          g2(n)%ni = ni *2-1; g2(n)%nj = nj*2-1;
          write(6,*) '  Size of block', n, 'ni =', ni, 'nj =', nj

          !allocate spatially varying dt, now inside the grid variable
          allocate(g(n)%dt(ni-1,nj-1))
          allocate(g2(n)%dt(g2(n)%ni-1,g2(n)%nj-1))
    
!         Allocate and then read the coordinates
          allocate(g(n)%x(ni,nj),g(n)%y(ni,nj))
          allocate(g2(n)%x(ni*2-1,nj*2-1),g2(n)%y(ni*2-1,nj*2-1))
          read(2) g(n)%x; read(2) g(n)%y;

!         Read the wall flag array
          allocate(g(n)%wall(ni,nj))
          allocate(g2(n)%wall(ni*2-1,nj*2-1))
          read(2) g(n)%wall;


          !Now actually generate the g2 mesh
          !placing the existing values in spread out indexes
            do i = 1, ni
                  do j = 1, nj 
                        g2(n)%x(2*i-1,2*j-1) = g(n)%x(i,j)
                        g2(n)%y(2*i-1,2*j-1) = g(n)%y(i,j)
                  end do
            end do  
            !filling in the new values at the old cell edges
            do i = 1, ni
                  do j = 1, nj-1 
                        g2(n)%x(2*i-1,2*j) = (g(n)%x(i,j)+g(n)%x(i,j+1))/2
                        g2(n)%y(2*i-1,2*j) = (g(n)%y(i,j)+g(n)%y(i,j+1))/2

                  end do
            end do  
            do i = 1, ni-1
                  do j = 1, nj

                        g2(n)%x(2*i,2*j-1) = (g(n)%x(i,j)+g(n)%x(i+1,j))/2
                        g2(n)%y(2*i,2*j-1) = (g(n)%y(i,j)+g(n)%y(i+1,j))/2

                  end do
            end do  
            !filling in the values at the cell centres
            do i = 1, ni-1
                  do j = 1, nj-1

                        g2(n)%x(2*i,2*j) = (g(n)%x(i,j)+g(n)%x(i+1,j)+g(n)%x(i,j+1)+g(n)%x(i+1,j+1))/4
                        g2(n)%y(2*i,2*j) = (g(n)%y(i,j)+g(n)%y(i+1,j)+g(n)%y(i,j+1)+g(n)%y(i+1,j+1))/4

                  end do
            end do  
            g2(n)%wall = .false.
            g2(n)%wall(:,[1,g2(n)%nj]) = .true.

!         Allocate all of the areas, side lengths and primary and secondary flow
!         variables and cell increments that are needed later
          allocate(g(n)%area(ni-1,nj-1),g(n)%lx_i(ni,nj-1),g(n)%ly_i(ni,nj-1), &
              g(n)%lx_j(ni-1,nj),g(n)%ly_j(ni-1,nj),g(n)%l_min(ni-1,nj-1))
          allocate(g2(n)%area(g2(n)%ni-1,g2(n)%nj-1),g2(n)%lx_i(g2(n)%ni,g2(n)%nj-1),g2(n)%ly_i(g2(n)%ni,g2(n)%nj-1), &
          g2(n)%lx_j(g2(n)%ni-1,g2(n)%nj),g2(n)%ly_j(g2(n)%ni-1,g2(n)%nj), g2(n)%l_min(g2(n)%ni-1,g2(n)%nj-1))

          allocate(g(n)%ro(ni,nj),g(n)%rovx(ni,nj),g(n)%rovy(ni,nj), &
              g(n)%roe(ni,nj))
          allocate(g2(n)%ro(ni*2-1,nj*2-1),g2(n)%rovx(ni*2-1,nj*2-1),g2(n)%rovy(ni*2-1,nj*2-1),g2(n)%roe(ni*2-1,nj*2-1))

          allocate(g(n)%dro(ni-1,nj-1),g(n)%drovx(ni-1,nj-1), &
              g(n)%drovy(ni-1,nj-1),g(n)%droe(ni-1,nj-1))
          allocate(g2(n)%dro(g2(n)%ni-1,g2(n)%nj-1),g2(n)%drovx(g2(n)%ni-1,g2(n)%nj-1), &
          g2(n)%drovy(g2(n)%ni-1,g2(n)%nj-1),g2(n)%droe(g2(n)%ni-1,g2(n)%nj-1))

          allocate(g(n)%p(ni,nj),g(n)%hstag(ni,nj),g(n)%vx(ni,nj), &
              g(n)%vy(ni,nj))
          allocate(g2(n)%p(ni*2-1,nj*2-1),g2(n)%hstag(ni*2-1,nj*2-1),g2(n)%vx(ni*2-1,nj*2-1),g2(n)%vy(ni*2-1,nj*2-1))

      end do

!     Read which blocks the inlet and outlets are on
      read(2) bcs%n_in, bcs%n_out
      write(6,*) '  Inlet on block', bcs%n_in, 'outlet on block', bcs%n_out

!     Allocate arrays to store static conditions at the inlet plane now we know
!     which block the inlet plane is on and its size in the j-direction
      allocate(bcs%ro(g(bcs%n_in)%nj),bcs%p(g(bcs%n_in)%nj))

!     Read the number of matching patches
      read(2) av%nm; allocate(p(av%nm));
      write(6,*) '  Number of matching patches nm =', av%nm

!     Loop over all matching patches and read the indices
      do m = 1,av%nm 

!         Read the length of the patch and the block indices
          read(2) p(m)%nk, p(m)%n_1, p(m)%n_2 
          write(6,*) '  Size of patch', m, 'nk =', p(m)%nk, 'joins block', &
              p(m)%n_1, 'and', p(m)%n_2

!         Allocate the length of the index lists
          allocate(p(m)%i_1(p(m)%nk),p(m)%j_1(p(m)%nk),p(m)%i_2(p(m)%nk), &
              p(m)%j_2(p(m)%nk))

!         Read the index lists
          read(2) p(m)%i_1,p(m)%j_1,p(m)%i_2,p(m)%j_2

      end do

!     Close the unit now everything has been read
      close(2)
      write(6,*)

      end subroutine read_mesh

end module


