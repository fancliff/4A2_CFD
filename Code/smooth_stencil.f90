
      module smooth_stencil

!     Packaging a subroutine in a module allows it to recieve the data
!     conveniently as assumed shape arrays
      
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine smooth_array(av,prop,prop_ref)

!     This subroutine smooths "prop" to stabilise the calculation, the basic 
!     solver uses second order smoothing, many improvements are possible.

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(inout) :: prop(:,:)
      real, intent(in) :: prop_ref
      real, dimension(size(prop,1),size(prop,2)) :: prop_avg_2, prop_avg_4_i, prop_avg_4_j, prop_avg_4, sfac_loc
      integer :: ni, nj
      real :: sf2, sf4
      logical :: local_smooth, fourth_smooth
      
!     sf2 and sf4 are only used if 4th_smooth is true. 
!     sf2 and sf4 must add up to 1
      sf2 = 0.2
      sf4 = 0.8
      fourth_smooth = .false.
      local_smooth = .false.

!     Get the block size and store locally for convenience
      ni = size(prop,1); nj = size(prop,2)

!     2ND ORDER SMOOTHING

!     Calculate the average values at the nodes in the interior region of the
!     mesh, use the four neighbouring nodes in the plus and minus i and 
!     j-directions.
      prop_avg_2(2:ni-1,2:nj-1) = ( prop(3:ni,2:nj-1) + prop(1:ni-2,2:nj-1) &
                                + prop(2:ni-1,3:nj) + prop(2:ni-1,1:nj-2) ) / 4.0

!     Edge values are also averaged in both the i and j-directions. Parallel to
!     the boundary the averaging is centred, the averages of two nodes are taken
!     either side of the current point. Perpendicular to the boundary the
!     algorithm is one-sided, the value at the current point is extrapolated
!     from the values at two nodes away from the boundary point.
      prop_avg_2(1,2:nj-1) = ( prop(1,1:nj-2) + prop(1,3:nj) &
                           + 2.0*prop(2,2:nj-1) - prop(3,2:nj-1) ) / 3.0

      prop_avg_2(ni,2:nj-1) = ( prop(ni,1:nj-2) + prop(ni,3:nj) &
                           + 2.0*prop(ni-1,2:nj-1) - prop(ni-2,2:nj-1) ) / 3.0                   

      prop_avg_2(2:ni-1,1) = ( prop(3:ni,1) + prop(1:ni-2,1) &
                           + 2.0*prop(2:ni-1,2) - prop(2:ni-1,3) ) / 3.0

      prop_avg_2(2:ni-1,nj) = ( prop(3:ni,nj) + prop(1:ni-2,nj) &
                           + 2.0*prop(2:ni-1,nj-1) - prop(2:ni-1,nj-2) ) / 3.0

!     The corner values are not currently smoothed
      prop_avg_2([1,ni],[1,nj]) = prop([1,ni],[1,nj])


      
!     4TH ORDER SMOOTHING
!     Only calculate if needed by checking fourth_smooth
!     Work out i direction first then add j direction
      
      if (fourth_smooth) then

!     i direction - central:
          prop_avg_4_i(3:ni-2,1:nj) = - 1/6 * prop(1:ni-4,1:nj) & 
                                      + 2/3 * prop(2:ni-3,1:nj) &
                                      + 2/3 * prop(4:ni-1,1:nj) &
                                      - 1/6 * prop(5:ni  ,1:nj)

!     j direction - central:
          prop_avg_4_j(1:ni,3:nj-2) = - 1/6 * prop(1:ni,1:nj-4) &
                                      + 2/3 * prop(1:ni,2:nj-3) &
                                      + 2/3 * prop(1:ni,4:nj-1) &
                                      - 1/6 * prop(1:ni,5:nj  )

!      i direction - semi one-sided:
          prop_avg_4_i(2,1:nj) = + 1/4 * prop(1,1:nj) &
                                 + 3/2 * prop(3,1:nj) &
                                 - 1   * prop(4,1:nj) &
                                 + 1/4 * prop(5,1:nj)

          prop_avg_4_i(ni-1,1:nj) = + 1/4 * prop(ni,1:nj)   &
                                    + 3/2 * prop(ni-2,1:nj) &
                                    - 1   * prop(ni-3,1:nj) &
                                    + 1/4 * prop(ni-4,1:nj)                                 

!      j direction - semi one-sided:
          prop_avg_4_j(1:ni,2) = + 1/4 * prop(1:ni,1) &
                                 + 3/2 * prop(1:ni,3) &
                                 - 1   * prop(1:ni,4) &
                                 + 1/4 * prop(1:ni,5)

          prop_avg_4_j(1:ni,nj-1) = + 1/4 * prop(1:ni,nj)   &
                                    + 3/2 * prop(1:ni,nj-2) &
                                    - 1   * prop(1:ni,nj-3) &
                                    + 1/4 * prop(1:ni,nj-4)

!     i direction - full one-sided (not corners):                                           
          prop_avg_4_i(1,2:nj-1) = + 4 * prop(2,2:nj-1) &
                                   - 6 * prop(3,2:nj-1) &
                                   + 4 * prop(4,2:nj-1) &
                                   - 1 * prop(5,2:nj-1)
                                   
          prop_avg_4_i(ni,2:nj-1) = + 4 * prop(ni-1,2:nj-1) &
                                    - 6 * prop(ni-2,2:nj-1) &
                                    + 4 * prop(ni-3,2:nj-1) &
                                    - 1 * prop(ni-4,2:nj-1)
                                    
!     j direction - full one-sided (not corners):                                           
          prop_avg_4_j(2:ni-1,1) = + 4 * prop(2:ni-1,2) &
                                   - 6 * prop(2:ni-1,3) &
                                   + 4 * prop(2:ni-1,4) &
                                   - 1 * prop(2:ni-1,5)
                                   
          prop_avg_4_j(2:ni-1,nj) = + 4 * prop(2:ni-1,2) &
                                    - 6 * prop(2:ni-1,3) &
                                    + 4 * prop(2:ni-1,4) &
                                    - 1 * prop(2:ni-1,5)

!     Add the i and j direction contributions together
          prop_avg_4 = prop_avg_4_i + prop_avg_4_j
      
!     Corners currently not smoothed
          prop_avg_4([1,ni],[1,nj]) = prop([1,ni],[1,nj]) 
          
      end if



!     Now apply the artificial viscosity by smoothing "prop" towards "prop_avg",
!     take (1-sfac) * the calculated value of the property + sfac * the average 
!     of the surrounding values. 
      
      if (local_smooth) then
          sfac_loc = av%sfac * abs(prop - prop_avg_2) / prop_ref
      else
          sfac_loc(:,:) = av%sfac
      end if
      
      if (fourth_smooth) then
          prop = (1 - sfac_loc) * prop + sfac_loc * &
                                ( sf2 * prop_avg_2 + sf4 * prop_avg_4 )
      else
          prop = (1 - sfac_loc) * prop + sfac_loc * prop_avg_2
      end if

      end subroutine smooth_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module smooth_stencil


