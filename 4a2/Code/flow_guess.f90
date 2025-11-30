      
module guess
      use types
      implicit none
contains
      subroutine flow_guess(av,g,bcs,guesstype)

!     This calculates an initial guess of the primary flowfield variables and
!     stores them at the nodes within the mesh dataytype

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(inout) :: g(:)
      type(t_bconds), intent(in) :: bcs
      integer, intent(in) :: guesstype
      integer :: i, j, ni, nj, j_mid, n
      
!     Variables required for the crude guess
      real :: t_out, v_out, ro_out, lx, ly, l


!     Variables required for the improved guess, you will need to add to these
!     INSERTED
      real :: mdot_out, mdot_block
      real :: t_lim
      real :: mach_lim
      real :: li_out(g(bcs%n_out)%ni)
      real, allocatable :: v_guess(:)
      real, allocatable :: ro_guess(:)
      real, allocatable :: t_guess(:)
      real, allocatable :: l_i(:)

      do n = 1, av%nn !probably need the outlet first and then go back through the domain
            allocate(v_guess(g(n)%ni))
            allocate(ro_guess(g(n)%ni))
            allocate(t_guess(g(n)%ni))
            allocate(l_i(g(n)%ni))
      !     Get the size of the mesh and store locally for convenience
            ni = g(n)%ni; nj = g(n)%nj;

      !     Assuming isentropic flow to the the exit plane calculate the static
      !     temperature and the exit velocity
            t_out = bcs%tstag * (bcs%p_out / bcs%pstag)**av%fgam
            v_out = (2 * av%cp * (bcs%tstag - t_out))**0.5
            ro_out = bcs%p_out / (av%rgas * t_out)

      !     Determine which guess calcation method to use by the value of "guesstype"
            if(guesstype == 1) then

            !         Store the exit density and internal energy as if they were uniform 
                  g(n)%ro = ro_out 
                  g(n)%roe  = ro_out * (av%cv * t_out + 0.5 * v_out**2)

            !         Calculate the gradient of the mesh lines in the centre of the domain
            !         to determine the assumed direction of the flow
                  j_mid = nj / 2
                  do i = 1,ni-1
                        lx = g(n)%lx_j(i,j_mid); ly = g(n)%ly_j(i,j_mid); 
                        l = hypot(lx,ly)
                        g(n)%rovx(i,:) = g(n)%ro(i,:) * v_out * abs(ly / l)  !added abs to stop flow reversal on the aerofoil block
                        g(n)%rovy(i,:) = -g(n)%ro(i,:) * v_out * lx / l
                  end do

            !         Copy the values to the "i = ni" nodes as an approximation
                  g(n)%rovx(ni,:) = g(n)%rovx(ni-1,:)
                  g(n)%rovy(ni,:) = g(n)%rovy(ni-1,:)

            !         Print the guess that has been calculated
                  write(6,*) 'Crude flow guess calculated'
                  write(6,*) '  At first point ro =', g(n)%ro(1,1), 'roe =', &
                        g(n)%roe(1,1), 'rovx =', g(n)%rovx(1,1), 'rovy =', g(n)%rovy(1,1)
                  write(6,*)

            else if(guesstype == 2) then 

      !         Calculate the length of each "i = const" line between the "j = 1" and 
      !         "j = nj" boundaries of the domain and store it in the local variable
      !         "l_i". You could calculate the length along each i-facet from the x 
      !         and y projected lengths with "hypot" and then sum them up in the
      !         second dimension with "sum". 
      !         INSERTED
                  l_i = sum(hypot(g(n)%ly_i, g(n)%lx_i),2)

      !         Use the exit temperature, density and velocity calculated for the 
      !         crude guess with "l_i" to estimate the mass flow rate at the exit
      !         INSERTED
                  li_out = sum(hypot(g(bcs%n_out)%ly_i, g(bcs%n_out)%lx_i),2)
                  mdot_out = ro_out * v_out * li_out(g(bcs%n_out)%ni) !!change this to either be exit mass flow for end block
                  !or calculated from continuity
                  mdot_block = mdot_out * l_i(ni)/li_out(g(bcs%n_out)%ni)

      !         Set a limit to the maximum allowable mach number in the initial
      !         guess, call this "mach_lim", calculate the corresponding temperature,
      !         called "t_lim"
      !         INSERTED
                  mach_lim = 0.5 !1
                  t_lim = bcs%tstag / (1+((av%gam-1)/2)*mach_lim**2)


      !         Now estimate the velocity and density at every "i = const" line, call 
      !         the velocity "v_guess(i)" and the density "ro_guess(i)":
      !             1. Assume density is constant at the exit value
      !             2. Use continuity and "l_i(i)" to estimate the velocity
      !             3. Assume stagnation temperature is constant for static temp
      !             4. Limit the static temperature, lookup intrinsic "max"
      !             5. Calculate the density throughout "ro_guess(i)"
      !             6. Update the estimate of the velocity "v_guess(i)" 
      !         INSERTED
                  ro_guess = ro_out
                  v_guess = mdot_block/(ro_guess*l_i) !l_i needs to be for all blocks
                  t_guess = bcs%tstag- (v_guess**2)/(2*av%cp)
                  t_guess = max(t_guess, t_lim)
                  ro_guess = bcs%rostag * (1+((av%gam-1)/2)*(v_guess**2/(av%gam*av%rgas*t_guess))) ** (-1/(av%gam-1))
                  v_guess = mdot_block/(ro_guess*l_i)



      !         Direct the calculated velocity to be parallel to the "j = const"
      !         gridlines for all values of i and j. This can be achieved with a 
      !         similar calculation to the "j = nj/2" one that was performed in the 
      !         crude guess. Then set all of ro, roe, rovx and rovy, note that roe 
      !         includes the kinetic energy component of the internal energy.
      !         INSERTED
            
                  do i = 1,ni-1 
                        do j = 1,nj
                              lx = g(n)%lx_j(i,j); ly = g(n)%ly_j(i,j); 
                              l = hypot(lx,ly)

                              g(n)%ro(i, :) = ro_guess(i)
                              g(n)%roe(i,j) =  g(n)%ro(i,j) * ( av%cv*t_guess(i) + 0.5*v_guess(i)**2 ) !t_out?
                              g(n)%rovx(i,j) = g(n)%ro(i,j) * abs(v_guess(i) * ly/l )
                              g(n)%rovy(i,j) = g(n)%ro(i,j) * (-v_guess(i) * lx/l)
                        end do
                  end do

            
      !         Make sure the guess has been copied for the "i = ni" values too
      !         INSERTED
                  g(n)%ro(ni,:) = g(n)%ro(ni-1,:)
                  g(n)%roe(ni,:) = g(n)%roe(ni-1,:)
                  g(n)%rovx(ni,:) = g(n)%rovx(ni-1,:)
                  g(n)%rovy(ni,:) = g(n)%rovy(ni-1,:)

      !         Print the first elements of the guess like for the crude guess
      !         INSERTED
                  write(6,*) 'Improved flow guess calculated'
                  write(6,*) '  At first point ro =', g(n)%ro(1,1), 'roe =', &
                  g(n)%roe(1,1), 'rovx =', g(n)%rovx(1,1), 'rovy =', g(n)%rovy(1,1)
                  write(6,*)

            end if

      !     The initial guess values derived from the boundary conditions are also
      !     useful as a reference to non-dimensionalise the convergence
            av%ro_ref = sum(g(n)%ro(1,:)) / nj
            av%roe_ref = sum(g(n)%roe(1,:)) / nj
            av%rov_ref = max(sum(g(n)%rovx(1,:)),sum(g(n)%rovy(1,:))) / nj

            deallocate(v_guess, ro_guess, t_guess, l_i)
            write(*,*) allocated(g(n)%rovx)
      end do

      end subroutine flow_guess 

end module


