      
      subroutine read_settings(av,bcs)

!     Read in the application variables and gas constants, they are in the
!     standard input file which is already assigned to unit 5

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(out) :: av
      type(t_bconds), intent(out) :: bcs
      character(len=80) :: tempname
      integer :: num_lines, i
      integer :: ios

!     Read the case name and trim to the required length
      read(5,*) tempname
      av%casename = trim(tempname)

!     You should read in the following variables sequentially and store them in
!     the dervived "av" datatype with the % syntax:
!         gam, rgas
!         cfl, sfac, d_max, facsec (for Crocco method), nrkuts (for Runge-Kutta)
!         nsteps
!         ni, nj
      read(5,*) av%gam, av%rgas
      read(5,*) av%cfl, av%sfac, av%d_max, av%facsec, av%nrkuts, av%spatial_acc
      read(5,*) av%nsteps
      read(5,*) av%ni, av%nj

!     Calculate other gas constants used throughout the calculation
      av%cp = av%rgas * av%gam / (av%gam - 1.0)
      av%cv = av%cp / av%gam
      av%fgam = (av%gam - 1.0) / av%gam

!     Scale the smoothing factor and convergence limits by cfl number, the 
!     changes over a timestep should be proportional to dt, which is 
!     proportional to cfl
      av%sfac = av%sfac * av%cfl
      av%d_max = av%d_max * av%cfl

!     Average convergence limit on residuals is set to half of the maximum
      av%d_avg = 0.5 * av%d_max

!     Read the inlet boundary condition data and store into the "bcs" datatype
!         pstag, tstag, alpha, rfin
      read(5,*) bcs%pstag, bcs%tstag, bcs%alpha, bcs%rfin      

!     Convert the inlet angle to radians
      bcs%alpha = bcs%alpha * 3.14159 / 180.0

!     Read the outlet static pressure and store into the "bcs" datatype
      read(5,*) bcs%p_out
      
!     Time varying-boundary condition stuff
!     Open the time varying boundary conditions and assign to unit 10

      ! First pass to determine number of lines in the file
      num_lines = 0
      open(unit=10, file='bcs_' // av%casename // '.txt', status="old", action="read", iostat=ios)
      if (ios /= 0) then
          print *, "Error opening file: ", 'bcs_' // av%casename // '.txt'
          stop
      endif

      do
          read(10, *, iostat=ios)
          if (ios /= 0) exit
          num_lines = num_lines + 1
      end do
      
      ! Rewind the file pointer to the beginning
      rewind(10)

      ! Allocate arrays based on the number of lines
      allocate(bcs%t_in(num_lines), bcs%pstag_in(num_lines), bcs%p_out_in(num_lines))

      ! Read the file again to fill the arrays
      do i = 1, num_lines
          read(10, *) bcs%t_in(i), bcs%pstag_in(i), bcs%p_out_in(i)
      end do
      close(10)


!     Print the settings to check they have been read, you can use this syntax
!     anywhere else you want in the program to debug your code
      write(6,*)
      write(6,*) 'Solver begins on ', av%casename, ' case'
      write(6,*)
      write(6,*) 'Read application variables from file'
      write(6,*) '  rgas =', av%rgas, 'cp =', av%cp, 'cv =', av%cv
      write(6,*) '  CFL =', av%cfl, 'sfac =', av%sfac
      write(6,*) '  Convergence  d_max =', av%d_max
      write(6,*) '  facsec =', av%facsec, 'nrkuts =', av%nrkuts, & 
          'spatial accuracy =', av%spatial_acc
      write(6,*) '  Mesh size  ni =', av%ni, 'nj =', av%nj
      write(6,*) '  Inlet  pstag =', bcs%pstag, 'tstag =', bcs%tstag, &
          'alpha =', bcs%alpha
      write(6,*) '  Outlet  p_out =', bcs%p_out
      write(6,*)

      end subroutine read_settings


