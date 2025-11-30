
module euler 
      use types
      use flux_stencil
      use smooth_stencil
      implicit none
contains
      subroutine euler_iteration(av,g)

!     This subroutine calculates the fluxes into each cell and then sums them to
!     update the primary flow properties

!     Explicitly declare the required variables
      use types
      use flux_stencil
      use smooth_stencil
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g(:)
      real, allocatable :: mass_i(:,:), flux_i(:,:)
      real, allocatable :: mass_j(:,:), flux_j(:,:)
      real, allocatable :: avg_ihstag(:,:), avg_ip(:,:), avg_irovx(:,:), avg_irovy(:,:), avg_ivx(:,:), avg_ivy(:,:)
      real, allocatable :: avg_jhstag(:,:), avg_jp(:,:), avg_jrovx(:,:), avg_jrovy(:,:), avg_jvx(:,:), avg_jvy(:,:)

      integer :: i, j, ni, nj, n 

      do n = 1,av%nn
            
            allocate(mass_i(g(n)%ni,g(n)%nj-1), flux_i(g(n)%ni,g(n)%nj-1))
            allocate(mass_j(g(n)%ni-1,g(n)%nj), flux_j(g(n)%ni-1,g(n)%nj))
            allocate(avg_ihstag(g(n)%ni,g(n)%nj-1), avg_ip(g(n)%ni,g(n)%nj-1), avg_irovx(g(n)%ni,g(n)%nj-1)&
            &, avg_irovy(g(n)%ni,g(n)%nj-1), avg_ivx(g(n)%ni,g(n)%nj-1), avg_ivy(g(n)%ni,g(n)%nj-1))
            allocate(avg_jhstag(g(n)%ni-1,g(n)%nj), avg_jp(g(n)%ni-1,g(n)%nj), avg_jrovx(g(n)%ni-1,g(n)%nj),&
            & avg_jrovy(g(n)%ni-1,g(n)%nj), avg_jvx(g(n)%ni-1,g(n)%nj), avg_jvy(g(n)%ni-1,g(n)%nj))

      !     Get the block size and store locally for convenience
            ni = g(n)%ni; nj = g(n)%nj

      !     Setup the continuity equation by calculating the mass flow through
      !     the facets in both the i and j-directions. Store these values in
      !     "mass_i" and "mass_j"
      !     INSERTED
            avg_jrovx =  0.5 * (g(n)%rovx(1:ni-1,:)+g(n)%rovx(2:ni,:))
            avg_jrovy =  0.5 * (g(n)%rovy(1:ni-1,:)+g(n)%rovy(2:ni,:))
            avg_irovx =  0.5 * (g(n)%rovx(:,1:nj-1)+g(n)%rovx(:,2:nj))
            avg_irovy =  0.5 * (g(n)%rovy(:,1:nj-1)+g(n)%rovy(:,2:nj))
            mass_i = avg_irovx*g(n)%lx_i + avg_irovy*g(n)%ly_i
            mass_j = avg_jrovx*g(n)%lx_j + avg_jrovy*g(n)%ly_j
      
      !     Apply the wall boundary condition by checking that two nodes at the
      !     end of a facet are both on a wall, if so then the appropriate mass
      !     flow array is set to have zero flow through that facet
            where(g(n)%wall(1:ni-1,:) .and. g(n)%wall(2:ni,:)) mass_j = 0 
            where(g(n)%wall(:,1:nj-1) .and. g(n)%wall(:,2:nj)) mass_i = 0 

      !     Update the density with mass fluxes by calling "sum_fluxes"
      !     INSERTED
            call sum_fluxes(av, mass_i, mass_j, g(n)%area, g(n)%ro, g(n)%dro,g(n)%dt)

      !     Setup the conservation of energy equation by calculated the enthalpy flux
      !     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
      !     and "mass_j" from before
      !     INSERTED
            avg_jhstag =  0.5 * (g(n)%hstag(1:ni-1,:)+g(n)%hstag(2:ni,:))
            avg_ihstag =  0.5 * (g(n)%hstag(:,1:nj-1)+g(n)%hstag(:,2:nj))
            flux_i = avg_ihstag * mass_i
            flux_j = avg_jhstag * mass_j

      !     Update the internal energy with enthalpy fluxes
      !     INSERTED
            call sum_fluxes(av, flux_i, flux_j, g(n)%area, g(n)%roe, g(n)%droe,g(n)%dt)

      !     Setup the x-momentum equation including momentum flux and pressure forces
      !     INSERTED
            avg_jp =  0.5 * (g(n)%p(1:ni-1,:)+g(n)%p(2:ni,:))
            avg_ip =  0.5 * (g(n)%p(:,1:nj-1)+g(n)%p(:,2:nj))
            avg_jvx =  0.5 * (g(n)%vx(1:ni-1,:)+g(n)%vx(2:ni,:))
            avg_ivx =  0.5 * (g(n)%vx(:,1:nj-1)+g(n)%vx(:,2:nj))
            flux_i = avg_ivx * mass_i + avg_ip * g(n)%lx_i
            flux_j = avg_jvx * mass_j + avg_jp * g(n)%lx_j

      !     Update the x-momentum with momentum flux
      !     INSERTED
            call sum_fluxes(av, flux_i, flux_j, g(n)%area, g(n)%rovx, g(n)%drovx,g(n)%dt)

      !     Setup the y-momentum equation including momentum flux and pressure forces
      !     INSERTED
            avg_jvy =  0.5 * (g(n)%vy(1:ni-1,:)+g(n)%vy(2:ni,:))
            avg_ivy =  0.5 * (g(n)%vy(:,1:nj-1)+g(n)%vy(:,2:nj))
            flux_i = avg_ivy * mass_i + avg_ip * g(n)%ly_i
            flux_j = avg_jvy * mass_j + avg_jp * g(n)%ly_j

      !     Update the y-momentum with momentum flux
      !     INSERTED
            call sum_fluxes(av, flux_i, flux_j, g(n)%area, g(n)%rovy, g(n)%drovy,g(n)%dt)

      !     Add artificial viscosity by smoothing all of the primary flow variables
            call smooth_array(av,g(n)%ro)
            call smooth_array(av,g(n)%roe)
            call smooth_array(av,g(n)%rovx)
            call smooth_array(av,g(n)%rovy)

            deallocate(mass_i, flux_i)
            deallocate( mass_j, flux_j)
            deallocate( avg_ihstag, avg_ip, avg_irovx, avg_irovy, avg_ivx, avg_ivy)
            deallocate( avg_jhstag, avg_jp, avg_jrovx, avg_jrovy, avg_jvx, avg_jvy)
      end do

      end subroutine euler_iteration
end module


