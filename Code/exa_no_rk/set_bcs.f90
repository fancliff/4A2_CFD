      
      subroutine set_bcs(av,bcs)

!     This subroutine sets up the time varying boundary condition arrays
!     By interpolating the input values from read_settings

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_bconds), intent(inout) :: bcs
      integer :: i
      
      ! Allocate actual bcs arrays based on number of steps
      allocate(bcs%t_var(av%nsteps), bcs%pstag_var(av%nsteps), bcs%p_out_var(av%nsteps))
      
      ! Generate time array to be used for interpolation of data
      do i=1,av%nsteps
          bcs%t_var(i) = (i-1) * av%dt
      end do
      
      ! Interpolate the input data to match the number of solver timesteps
      call interp(bcs%t_in,bcs%pstag_in,bcs%t_var,bcs%pstag_var)
      call interp(bcs%t_in,bcs%p_out_in,bcs%t_var,bcs%p_out_var)

      end subroutine set_bcs
