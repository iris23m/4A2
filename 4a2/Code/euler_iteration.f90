
      subroutine euler_iteration(av,g)

!     This subroutine calculates the fluxes into each cell and then sums them to
!     update the primary flow properties

!     Explicitly declare the required variables
      use types
      use flux_stencil
      use smooth_stencil
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      real, dimension(g%ni,g%nj-1) :: mass_i, flux_i
      real, dimension(g%ni-1,g%nj) :: mass_j, flux_j
      integer :: i, j, ni, nj

      real, dimension(g%ni,g%nj-1) :: avg_ihstag, avg_ip, avg_irovx, avg_irovy, avg_ivx, avg_ivy
      real, dimension(g%ni-1,g%nj) :: avg_jhstag, avg_jp, avg_jrovx, avg_jrovy, avg_jvx, avg_jvy

!     Get the block size and store locally for convenience
      ni = g%ni; nj = g%nj

!     Setup the continuity equation by calculating the mass flow through
!     the facets in both the i and j-directions. Store these values in
!     "mass_i" and "mass_j"
!     INSERTED
      avg_jrovx =  0.5 * (g%rovx(1:ni-1,:)+g%rovx(2:ni,:))
      avg_jrovy =  0.5 * (g%rovy(1:ni-1,:)+g%rovy(2:ni,:))
      avg_irovx =  0.5 * (g%rovx(:,1:nj-1)+g%rovx(:,2:nj))
      avg_irovy =  0.5 * (g%rovy(:,1:nj-1)+g%rovy(:,2:nj))
      mass_i = avg_irovx*g%lx_i + avg_irovy*g%ly_i
      mass_j = avg_jrovx*g%lx_j + avg_jrovy*g%ly_j
     
!     Apply the wall boundary condition by checking that two nodes at the
!     end of a facet are both on a wall, if so then the appropriate mass
!     flow array is set to have zero flow through that facet
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 

!     Update the density with mass fluxes by calling "sum_fluxes"
!     INSERTED
      call sum_fluxes(av, mass_i, mass_j, g%area, g%ro, g%dro)

!     Setup the conservation of energy equation by calculated the enthalpy flux
!     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!     and "mass_j" from before
!     INSERTED
      avg_jhstag =  0.5 * (g%hstag(1:ni-1,:)+g%hstag(2:ni,:))
      avg_ihstag =  0.5 * (g%hstag(:,1:nj-1)+g%hstag(:,2:nj))
      flux_i = avg_ihstag * mass_i
      flux_j = avg_jhstag * mass_j

!     Update the internal energy with enthalpy fluxes
!     INSERTED
      call sum_fluxes(av, flux_i, flux_j, g%area, g%roe, g%droe)

!     Setup the x-momentum equation including momentum flux and pressure forces
!     INSERTED
      avg_jp =  0.5 * (g%p(1:ni-1,:)+g%p(2:ni,:))
      avg_ip =  0.5 * (g%p(:,1:nj-1)+g%p(:,2:nj))
      avg_jvx =  0.5 * (g%vx(1:ni-1,:)+g%vx(2:ni,:))
      avg_ivx =  0.5 * (g%vx(:,1:nj-1)+g%vx(:,2:nj))
      flux_i = avg_ivx * mass_i + avg_ip * g%lx_i
      flux_j = avg_jvx * mass_j + avg_jp * g%lx_j

!     Update the x-momentum with momentum flux
!     INSERTED
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovx, g%drovx)

!     Setup the y-momentum equation including momentum flux and pressure forces
!     INSERTED
      avg_jvy =  0.5 * (g%vy(1:ni-1,:)+g%vy(2:ni,:))
      avg_ivy =  0.5 * (g%vy(:,1:nj-1)+g%vy(:,2:nj))
      flux_i = avg_ivy * mass_i + avg_ip * g%ly_i
      flux_j = avg_jvy * mass_j + avg_jp * g%ly_j

!     Update the y-momentum with momentum flux
!     INSERTED
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovy, g%drovy)

!     Add artificial viscosity by smoothing all of the primary flow variables
      call smooth_array(av,g%ro)
      call smooth_array(av,g%roe)
      call smooth_array(av,g%rovx)
      call smooth_array(av,g%rovy)
      

      end subroutine euler_iteration


