module interp 
      use types
      implicit none
contains
subroutine interpolate_mesh(g,g2, bcs,av)
      use types 
      implicit none
      type(t_bconds), intent(inout) :: bcs
      type(t_grid), intent(inout) :: g(:)
      type(t_grid), intent(inout) :: g2(:)
      type(t_appvars), intent(inout) :: av
      integer :: ni, nj, i, j,n
      real, allocatable :: tempro(:)

      do n =1,av%nn
            ni = g(n)%ni  
            nj = g(n)%nj

            if (n == bcs%n_in) then
                  ! resize the bconds
                  allocate(tempro(2*nj-1))
                  do j = 1,nj
                  tempro(2*j-1) = bcs%ro(j)
                  end do 
                  do j = 1,nj-1
                  tempro(2*j) = (bcs%ro(j)+bcs%ro(j+1))/2
                  end do
                  deallocate(bcs%ro)
                  allocate(bcs%ro(2*nj-1))
                  bcs%ro = tempro
            end if

            !placing the existing values in spread out indexes
            do i = 1, ni
                  do j = 1, nj 
                        g2(n)%ro(2*i-1,2*j-1) = g(n)%ro(i,j)
                        g2(n)%roe(2*i-1,2*j-1) = g(n)%roe(i,j)
                        g2(n)%rovx(2*i-1,2*j-1) = g(n)%rovx(i,j)
                        g2(n)%rovy(2*i-1,2*j-1) = g(n)%rovy(i,j)
                  
                  end do
            end do  
            !filling in the new values at the old cell edges
            do i = 1, ni
                  do j = 1, nj-1 
                        g2(n)%ro(2*i-1,2*j) = (g(n)%ro(i,j)+g(n)%ro(i,j+1))/2
                        g2(n)%roe(2*i-1,2*j) = (g(n)%roe(i,j)+g(n)%roe(i,j+1))/2
                        g2(n)%rovx(2*i-1,2*j) = (g(n)%rovx(i,j)+g(n)%rovx(i,j+1))/2
                        g2(n)%rovy(2*i-1,2*j) = (g(n)%rovy(i,j)+g(n)%rovy(i,j+1))/2
                        
                  end do
            end do  
            do i = 1, ni-1
                  do j = 1, nj

                        g2(n)%ro(2*i,2*j-1) = (g(n)%ro(i,j)+g(n)%ro(i+1,j))/2
                        g2(n)%roe(2*i,2*j-1) = (g(n)%roe(i,j)+g(n)%roe(i+1,j))/2
                        g2(n)%rovx(2*i,2*j-1) = (g(n)%rovx(i,j)+g(n)%rovx(i+1,j))/2
                        g2(n)%rovy(2*i,2*j-1) = (g(n)%rovy(i,j)+g(n)%rovy(i+1,j))/2

                  end do
            end do  
            !filling in the values at the cell centres
            do i = 1, ni-1
                  do j = 1, nj-1

                        g2(n)%ro(2*i,2*j) = (g(n)%ro(i,j)+g(n)%ro(i+1,j)+g(n)%ro(i,j+1)+g(n)%ro(i+1,j+1))/4
                        g2(n)%roe(2*i,2*j) = (g(n)%roe(i,j)+g(n)%roe(i+1,j)+g(n)%roe(i,j+1)+g(n)%roe(i+1,j+1))/4
                        g2(n)%rovx(2*i,2*j) = (g(n)%rovx(i,j)+g(n)%rovx(i+1,j)+g(n)%rovx(i,j+1)+g(n)%rovx(i+1,j+1))/4
                        g2(n)%rovy(2*i,2*j) = (g(n)%rovy(i,j)+g(n)%rovy(i+1,j)+g(n)%rovy(i,j+1)+g(n)%rovy(i+1,j+1))/4

                  end do
            end do  
      end do

end subroutine interpolate_mesh

end module

