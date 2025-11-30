subroutine interpolate_mesh(g,g2, bcs)
      use types 
      implicit none
      type(t_bconds), intent(inout) :: bcs
      type(t_grid), intent(inout) :: g
      type(t_grid), intent(inout) :: g2
      integer :: ni, nj, i, j
      real, allocatable :: tempro(:)

     
      ni = g%ni  
      nj = g%nj

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

      !placing the existing values in spread out indexes
      do i = 1, ni
            do j = 1, nj 
                  g2%ro(2*i-1,2*j-1) = g%ro(i,j)
                  g2%roe(2*i-1,2*j-1) = g%roe(i,j)
                  g2%rovx(2*i-1,2*j-1) = g%rovx(i,j)
                  g2%rovy(2*i-1,2*j-1) = g%rovy(i,j)
                
            end do
      end do  
      !filling in the new values at the old cell edges
      do i = 1, ni
            do j = 1, nj-1 
                  g2%ro(2*i-1,2*j) = (g%ro(i,j)+g%ro(i,j+1))/2
                  g2%roe(2*i-1,2*j) = (g%roe(i,j)+g%roe(i,j+1))/2
                  g2%rovx(2*i-1,2*j) = (g%rovx(i,j)+g%rovx(i,j+1))/2
                  g2%rovy(2*i-1,2*j) = (g%rovy(i,j)+g%rovy(i,j+1))/2
                  
            end do
      end do  
      do i = 1, ni-1
            do j = 1, nj

                  g2%ro(2*i,2*j-1) = (g%ro(i,j)+g%ro(i+1,j))/2
                  g2%roe(2*i,2*j-1) = (g%roe(i,j)+g%roe(i+1,j))/2
                  g2%rovx(2*i,2*j-1) = (g%rovx(i,j)+g%rovx(i+1,j))/2
                  g2%rovy(2*i,2*j-1) = (g%rovy(i,j)+g%rovy(i+1,j))/2

            end do
      end do  
      !filling in the values at the cell centres
      do i = 1, ni-1
            do j = 1, nj-1

                  g2%ro(2*i,2*j) = (g%ro(i,j)+g%ro(i+1,j)+g%ro(i,j+1)+g%ro(i+1,j+1))/4
                  g2%roe(2*i,2*j) = (g%roe(i,j)+g%roe(i+1,j)+g%roe(i,j+1)+g%roe(i+1,j+1))/4
                  g2%rovx(2*i,2*j) = (g%rovx(i,j)+g%rovx(i+1,j)+g%rovx(i,j+1)+g%rovx(i+1,j+1))/4
                  g2%rovy(2*i,2*j) = (g%rovy(i,j)+g%rovy(i+1,j)+g%rovy(i,j+1)+g%rovy(i+1,j+1))/4

            end do
      end do  

end subroutine interpolate_mesh

