
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
      real, dimension(g%ni,g%nj-1) :: mass_i, flux_i, flux_i_2
      real, dimension(g%ni-1,g%nj) :: mass_j, flux_j, flux_j_2
      integer :: i, j, ni, nj
      logical :: fourth_flux
      
!     set the parameter to switch between fourth order and second order flux
      fourth_flux = .true.
!      fourth_flux = .false.

!     Get the block size and store locally for convenience
      ni = g%ni; nj = g%nj

!     Setup the continuity equation by calculating the mass flow through
!     the facets in both the i and j-directions. Store these values in
!     "mass_i" and "mass_j"
      call get_flux_i(g%rovx,flux_i,fourth_flux)
      call get_flux_i(g%rovy,flux_i_2,fourth_flux)
      call get_flux_j(g%rovx,flux_j,fourth_flux)
      call get_flux_j(g%rovy,flux_j_2,fourth_flux) 
      
      mass_i = flux_i * g%lx_i + flux_i_2 * g%ly_i
      mass_j = flux_j * g%lx_j + flux_j_2 * g%ly_j    
     
!     Apply the wall boundary condition by checking that two nodes at the
!     end of a facet are both on a wall, if so then the appropriate mass
!     flow array is set to have zero flow through that facet
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0.0
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0.0

!     Update the density with mass fluxes by calling "sum_fluxes"
      call sum_fluxes(av,mass_i,mass_j,g%area,g%ro,g%ro_start,g%dro)

!     Setup the conservation of energy equation by calculated the enthalpy flux
!     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!     and "mass_j" from before
      call get_flux_i(g%hstag,flux_i,fourth_flux)
      call get_flux_j(g%hstag,flux_j,fourth_flux)
      flux_i = flux_i * mass_i
      flux_j = flux_j * mass_j

!     Update the internal energy with enthalpy fluxes
      call sum_fluxes(av,flux_i,flux_j,g%area,g%roe,g%roe_start,g%droe)

!     Setup the x-momentum equation including momentum flux and pressure forces
      call get_flux_i(g%vx,flux_i,fourth_flux)
      call get_flux_i(g%p,flux_i_2,fourth_flux)
      call get_flux_j(g%vx,flux_j,fourth_flux)
      call get_flux_j(g%p,flux_j_2,fourth_flux) 
      
      flux_i = flux_i * mass_i + flux_i_2 * g%lx_i
      flux_j = flux_j * mass_j + flux_j_2 * g%lx_j

!     Update the x-momentum with momentum flux
      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovx,g%rovx_start,g%drovx)

!     Setup the y-momentum equation including momentum flux and pressure forces
      call get_flux_i(g%vy,flux_i,fourth_flux)
      call get_flux_i(g%p,flux_i_2,fourth_flux)
      call get_flux_j(g%vy,flux_j,fourth_flux)
      call get_flux_j(g%p,flux_j_2,fourth_flux) 
      
      flux_i = flux_i * mass_i + flux_i_2 * g%ly_i
      flux_j = flux_j * mass_j + flux_j_2 * g%ly_j

!     Update the y-momentum with momentum flux
      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovy,g%rovy_start,g%drovy)
            

!     Add artificial viscosity by smoothing all of the primary flow variables
!     Include reference property value for use in determining sfac_loc
      call smooth_array(av, g%ro, av%ro_ref)
      call smooth_array(av, g%roe, av%roe_ref)
      call smooth_array(av, g%rovx, av%rov_ref)
      call smooth_array(av, g%rovy, av%rov_ref)
      

      end subroutine euler_iteration


