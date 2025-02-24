      
      subroutine flow_guess(av,g,bcs,guesstype)

!     This calculates an initial guess of the primary flowfield variables and
!     stores them at the nodes within the mesh dataytype

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(in) :: bcs
      integer, intent(in) :: guesstype
      integer :: i, j, ni, nj, j_mid, i_mid
      
!     Variables required for the crude guess
      real :: t_out, v_out, ro_out, lx, ly, l

!     Variables required for the improved guess, you will need to add to these
      real, dimension(g%ni) :: l_i, v_guess, ro_guess, t_guess
      real :: mdot, mach_lim, t_lim
      real, dimension(g%ni-1,g%nj) :: dx, dy, dl
      
!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Assuming isentropic flow to the the exit plane calculate the static
!     temperature and the exit velocity
      t_out = bcs%tstag * (bcs%p_out / bcs%pstag)**av%fgam
      v_out = sqrt(2.0 * av%cp * (bcs%tstag - t_out))
      ro_out = bcs%p_out / (av%rgas * t_out)

!     Special flow guess for sod shock tube
      if(guesstype == 0) then
          i_mid = ni/2 !fortran integer division automatically floored (towards zero for -ve)
          g%ro(1:i_mid,:) = 1.0000
          g%ro(i_mid+1:ni,:) = 0.12500
          g%rovx(:,:) = 0.00000
          g%rovy(:,:) = 0.00000
          g%roe(1:i_mid,:)    = 1.0000 / (av%gam-1)
          g%roe(i_mid+1:ni,:) = 0.1000/ (av%gam-1)
      end if

!     Determine which guess calcation method to use by the value of "guesstype"
      if(guesstype == 1) then

!         Store the exit density and internal energy as if they were uniform 
          g%ro = ro_out 
          g%roe  = g%ro * (av%cv * t_out + 0.5 * v_out**2)

!         Calculate the gradient of the mesh lines in the centre of the domain
!         to determine the assumed direction of the flow
          j_mid = nj / 2
          do i = 1,ni-1
              lx = g%lx_j(i,j_mid); ly = g%ly_j(i,j_mid); 
              l = hypot(lx,ly)
              g%rovx(i,:) = g%ro(i,:) * v_out * ly / l
              g%rovy(i,:) = -g%ro(i,:) * v_out * lx / l
          end do

!         Copy the values to the "i = ni" nodes as an approximation
          g%rovx(ni,:) = g%rovx(ni-1,:)
          g%rovy(ni,:) = g%rovy(ni-1,:)

!         Print the guess that has been calculated
          write(6,*) 'Crude flow guess calculated'
          write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
              g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
          write(6,*)

      else if(guesstype == 2) then 

!         Calculate the length of each "i = const" line between the "j = 1" and 
!         "j = nj" boundaries of the domain and store it in the local variable
!         "l_i". You could calculate the length along each i-facet from the x 
!         and y projected lengths with "hypot" and then sum them up in the
!         second dimension with "sum". 
          l_i = sum( hypot(g%lx_i,g%ly_i) , 2)

!         Use the exit temperature, density and velocity calculated for the 
!         crude guess with "l_i" to estimate the mass flow rate at the exit
          mdot = ro_out * v_out * l_i(ni)

!         Set a limit to the maximum allowable mach number in the initial
!         guess, call this "mach_lim", calculate the corresponding temperature,
!         called "t_lim"
          mach_lim = 0.9999
          t_lim = bcs%tstag / (1.0 + 0.5 * (av%gam-1.0) * mach_lim**2)

!         Now estimate the velocity and density at every "i = const" line, call 
!         the velocity "v_guess(i)" and the density "ro_guess(i)":
!             1. Assume density is constant at the exit value
!             2. Use continuity and "l_i(i)" to estimate the velocity
!             3. Assume stagnation temperature is constant for static temp
!             4. Limit the static temperature, lookup intrinsic "max"
!             5. Calculate the density throughout "ro_guess(i)"
!             6. Update the estimate of the velocity "v_guess(i)" 
          v_guess = mdot / ( ro_out * l_i )
          t_guess = max(t_lim, (bcs%tstag - (v_guess**2)/(2.0*av%cp)))
          ro_guess = bcs%pstag * (t_guess/bcs%tstag)**(1.0/av%fgam) & 
          					/ (av%rgas * t_guess)
          v_guess = mdot / ( ro_guess * l_i )

!         Direct the calculated velocity to be parallel to the "j = const"
!         gridlines for all values of i and j. This can be achieved with a 
!         similar calculation to the "j = nj/2" one that was performed in the 
!         crude guess. Then set all of ro, roe, rovx and rovy, note that roe 
!         includes the kinetic energy component of the internal energy.
          dx = g%lx_j(1:ni-1,:)
          dy = g%ly_j(1:ni-1,:)
          dl = hypot(dx,dy)
          do i = 1,ni-1
              g%ro(i,:) = ro_guess(i)
              g%roe(i,:) = ro_guess(i) * (av%cv * t_guess(i) + 0.5 * v_guess(i)**2.0)
              g%rovx(i,:) = g%ro(i,:) * v_guess(i) * dy(i,:) / dl(i,:)
              g%rovy(i,:) = -g%ro(i,:) * v_guess(i) * dx(i,:) / dl(i,:)
          end do					
          
              
!         Make sure the guess has been copied for the "i = ni" values too
      	  g%ro(ni,:) = g%ro(ni-1,:)
      	  g%roe(ni,:) = g%roe(ni-1,:)
          g%rovx(ni,:) = g%rovx(ni-1,:)
          g%rovy(ni,:) = g%rovy(ni-1,:)

!         Print the first elements of the guess like for the crude guess
          write(6,*) 'Improved flow guess calculated'
          write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
              g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
          write(6,*)

      end if

!     The initial guess values derived from the boundary conditions are also
!     useful as a reference to non-dimensionalise the convergence
      av%ro_ref = sum(g%ro(1,:)) / nj
      av%roe_ref = sum(g%roe(1,:)) / nj
      
!     need a new rov_ref for sod shock tube as velocity is 0 everywhere initially
!     use max flow speed at final solution (t=0.2) which can be found in sod.raw

      if(guesstype == 0) then
          av%rov_ref = 0.89095233
      else
          av%rov_ref = max(sum(g%rovx(1,:)),sum(g%rovy(1,:))) / nj
      end if
      
      
      end subroutine flow_guess


