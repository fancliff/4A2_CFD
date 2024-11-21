
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
      real, dimension(size(prop,1),size(prop,2)) :: prop_avg_2
      real, dimension(size(prop,1),size(prop,2)) :: prop_avg_4_i
      real, dimension(size(prop,1),size(prop,2)) :: prop_avg_4_j
      real, dimension(size(prop,1),size(prop,2)) :: prop_avg_4
      real, dimension(size(prop,1),size(prop,2)) :: sfac_loc
      integer :: ni, nj
      real :: sf2, sf4
      
!     sf2 and sf4 are only used if 4th_smooth is true. 
      sf4 = av%sfac
      sf2 = sf4/4.0

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
      
      !edges second order
      prop_avg_4(:,:) = prop_avg_2(:,:)

!     i direction - central:
      prop_avg_4_i(3:ni-2,3:nj-2) = - 1.0/6.0 * prop(1:ni-4,3:nj-2) &
                                  + 2.0/3.0 * prop(2:ni-3,3:nj-2) &
                                  + 2.0/3.0 * prop(4:ni-1,3:nj-2) &
                                  - 1.0/6.0 * prop(5:ni  ,3:nj-2)



!     j direction - central:
      prop_avg_4_j(3:ni-2,3:nj-2) = - 1.0/6.0 * prop(3:ni-2,1:nj-4) &
                                  + 2.0/3.0 * prop(3:ni-2,2:nj-3) &
                                  + 2.0/3.0 * prop(3:ni-2,4:nj-1) &
                                  - 1.0/6.0 * prop(3:ni-2,5:nj  )

  

!     Add the i and j direction contributions together and divide by 2
      prop_avg_4(3:ni-2,3:nj-2) = (prop_avg_4_i(3:ni-2,3:nj-2) + prop_avg_4_j(3:ni-2,3:nj-2)) / 2.0
          



!     Now apply the artificial viscosity by smoothing "prop" towards "prop_avg",
!     take (1-sfac) * the calculated value of the property + sfac * the average 
!     of the surrounding values. 
      

!     sf2 = av%sfac * sf2
!     sf4 = av%sfac * sf4	  
      prop = (1 - sf2 -sf4) * prop + sf2 * prop_avg_2 + sf4 * prop_avg_4
!      prop = (1-sfac_loc)*prop + sfac_loc*prop_avg_4         
          

      end subroutine smooth_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module smooth_stencil


