module patch
    use types
    implicit none
contains

subroutine patch_blocks(g,p,av)
    use types
    implicit none
    type(t_match), intent(inout) :: p(:)
    type(t_grid), intent(inout) :: g(:)
    type(t_appvars), intent(inout) :: av
    integer :: n, b1,b2, m
    real :: roval, rovxval, rovyval, roeval, vxval, vyval, pval, hstagval

    do n = 1,av%nm
        b1 = p(n)%n_1
        b2 = p(n)%n_2

       
        do m = 1,p(n)%nk

            roval = (g(b1)%ro(p(n)%i_1(m),p(n)%j_1(m)) + g(b2)%ro(p(n)%i_2(m),p(n)%j_2(m)))/2
            g(b1)%ro(p(n)%i_1(m),p(n)%j_1(m)) = roval
            g(b2)%ro(p(n)%i_2(m),p(n)%j_2(m)) = roval

            rovxval = (g(b1)%rovx(p(n)%i_1(m),p(n)%j_1(m)) + g(b2)%rovx(p(n)%i_2(m),p(n)%j_2(m)))/2
            g(b1)%rovx(p(n)%i_1(m),p(n)%j_1(m)) = rovxval
            g(b2)%rovx(p(n)%i_2(m),p(n)%j_2(m)) = rovxval

            rovyval = (g(b1)%rovy(p(n)%i_1(m),p(n)%j_1(m)) + g(b2)%rovy(p(n)%i_2(m),p(n)%j_2(m)))/2
            g(b1)%rovy(p(n)%i_1(m),p(n)%j_1(m)) = rovyval
            g(b2)%rovy(p(n)%i_2(m),p(n)%j_2(m)) = rovyval

            roeval = (g(b1)%roe(p(n)%i_1(m),p(n)%j_1(m)) + g(b2)%roe(p(n)%i_2(m),p(n)%j_2(m)))/2
            g(b1)%roe(p(n)%i_1(m),p(n)%j_1(m)) = roeval
            g(b2)%roe(p(n)%i_2(m),p(n)%j_2(m)) = roeval

    !       average secondary variables at patch
            vxval = (g(b1)%vx(p(n)%i_1(m),p(n)%j_1(m)) + g(b2)%vx(p(n)%i_2(m),p(n)%j_2(m)))/2
            g(b1)%vx(p(n)%i_1(m),p(n)%j_1(m)) = vxval
            g(b2)%vx(p(n)%i_2(m),p(n)%j_2(m)) = vxval

            vyval = (g(b1)%vy(p(n)%i_1(m),p(n)%j_1(m)) + g(b2)%vy(p(n)%i_2(m),p(n)%j_2(m)))/2
            g(b1)%vy(p(n)%i_1(m),p(n)%j_1(m)) = vyval
            g(b2)%vy(p(n)%i_2(m),p(n)%j_2(m)) = vyval

            pval = (g(b1)%p(p(n)%i_1(m),p(n)%j_1(m)) + g(b2)%p(p(n)%i_2(m),p(n)%j_2(m)))/2
            g(b1)%p(p(n)%i_1(m),p(n)%j_1(m)) = pval
            g(b2)%p(p(n)%i_2(m),p(n)%j_2(m)) = pval

            hstagval = (g(b1)%hstag(p(n)%i_1(m),p(n)%j_1(m)) + g(b2)%hstag(p(n)%i_2(m),p(n)%j_2(m)))/2
            g(b1)%hstag(p(n)%i_1(m),p(n)%j_1(m)) = hstagval
            g(b2)%hstag(p(n)%i_2(m),p(n)%j_2(m)) = hstagval
        end do
    end do
end subroutine patch_blocks

end module