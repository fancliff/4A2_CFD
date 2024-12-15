      
      subroutine set_timestep(av,g,bcs)

!     This subroutine sets a single value for the time step based on the 
!     stagnation speed of sound and the minimum length scale of any element

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(in) :: g
      type(t_bconds), intent(in) :: bcs
      real :: astag, v_max

!     Calculate the stagnation speed of sound from the inlet stagnation
!     temperature and gas constants
      astag = sqrt(av%gam * av%rgas * bcs%tstag)

!     Assume that the maximum flow speed is also equal to "astag". This will 
!     be pessimistic for subsonic flows but may be optimistic for supersonic 
!     flows. In the latter case the length of the time step as determined by 
!     may need to be reduced by improving this routine or varying the CFL number
      v_max = astag * 2.0

!     Calculate the timestep using the CFL number and store it in "av%dt"
!      av%dt = av%cfl * g%l_min / v_max

!     For Sod case the timestep is just a function of nsteps
      av%dt = 0.2 / av%nsteps

!     Print the calculated timestep and some intermediate values
!      write(6,*) 'a(stag): ', astag
!      write(6,*) 'V max: ', v_max
!      write(6,*) 'CFL no. ', av%cfl
!      write(6,*) 'L min: ', 6g%l_min
!      write(6,*)
      write(6,*) 'Runtime: 0.2s'
      write(6,*) 'Num steps: ', av%nsteps
      write(6,*) 'Timestep: ', av%dt
      write(6,*)

      end subroutine set_timestep


