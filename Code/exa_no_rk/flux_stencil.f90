
      module flux_stencil

!     Packaging a subroutine in a module allows it to recieve the data
!     conveniently as assumed shape arrays
      
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine sum_fluxes(av,flux_i,flux_j,area,prop,dcell)

!     This subroutine sums the fluxes into each cell, calculates the change in 
!     the cell property inside, distributes the change to the four nodes of the
!     cell and then adds it onto the flow property

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(in) :: flux_i(:,:), flux_j(:,:), area(:,:)
      real, intent(inout) :: prop(:,:)
      real, intent(inout) :: dcell(:,:)
      real, dimension(size(prop,1),size(prop,2)) :: dnode
      real, dimension(size(dcell,1),size(dcell,2)) :: dcell_temp, dcell_first
      integer :: ni, nj

!     Get the block size and store locally for convenience
      ni = size(prop,1); nj = size(prop,2)

!     Store previous value of dcell to dcell_temp for Crocco method
      dcell_temp = dcell

!     Use the finite volume method to find the change in the variables "prop"
!     over the timestep "dt", save it in the array "dcell"
      dcell = ( av%dt / area ) * & 
         (flux_i(1:ni-1,:) - flux_i(2:ni,:) + flux_j(:,1:nj-1) - flux_j(:,2:nj))

         
!     Crocco Method
!     Retain first order time accurate residual for next time-step
      dcell_first = dcell
      dcell = (1 + av%facsec) * dcell - av%facsec * dcell_temp
      

!     Now distribute the changes equally to the four corners of each cell. Each 
!     interior grid point receives one quarter of the change from each of the 
!     four cells adjacent to it.
      dnode(2:ni-1,2:nj-1) = ( dcell(2:ni-1,2:nj-1) + dcell(1:ni-2,2:nj-1) & 
                             + dcell(2:ni-1,1:nj-2) + dcell(1:ni-2,1:nj-2) ) / 4.0

!     Bounding edge nodes do not have four adjacent cells and so must be treated
!     differently, they only recieve half the change from each of the two
!     adjacent cells. Distribute the changes for the "i = 1 & ni" edges as well
!     as the "j = 1 & nj" edges. 
      dnode(1,2:nj-1) = (dcell(1,2:nj-1) + dcell(1,1:nj-2)) / 2.0
      dnode(ni,2:nj-1) = (dcell(ni-1,2:nj-1) + dcell(ni-1,1:nj-2)) / 2.0
      
      dnode(2:ni-1,1) = (dcell(2:ni-1,1) + dcell(1:ni-2,1)) / 2.0
      dnode(2:ni-1,nj) = (dcell(2:ni-1,nj-1) + dcell(1:ni-2,nj-1)) / 2.0

!     Finally distribute the changes to be to the four bounding corner points, 
!     these receive the full change from the single cell of which they form one 
!     corner.
      dnode([1,ni],[1,nj]) = dcell([1,ni-1],[1,nj-1])

!     Update the solution by adding the changes at the nodes "dnode" to the flow
!     property "prop"
      prop = prop + dnode

!     Return dcell to first order accurate after use to calculate nodes
      dcell = dcell_first
      
      
      end subroutine sum_fluxes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      subroutine get_flux_i(av,prop,flux_i)
      
!     gets the i flux component using either fourth or second order estimate

      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(in) :: prop(:,:)
      real, intent(inout) :: flux_i(:,:)
      integer :: ni, nj
      
      ni = size(prop,1); nj = size(prop,2)
      
      if(av%spatial_acc == 6) then
          
          flux_i(:,3:nj-3) = ( 11.0 * prop(:,1:nj-5) - 93.0 * prop(:,2:nj-4) & 
                             + 802.0 * ( prop(:,3:nj-3) + prop(:,4:nj-2) ) &
                             - 93.0 * prop(:,5:nj-1) + 11.0 * prop(:,6:nj) ) / 1440.0
          flux_i(:,2) = ( - 27.0 * prop(:,1) + 637.0 * prop(:,2) & 
                        + 1022.0 * prop(:,3) - 258.0 * prop(:,4) &
                        + 77.0 * prop(:,5) - 11.0 * prop(:,6) ) / 1440.0
          flux_i(:,nj-2) = ( - 27.0 * prop(:,nj) + 637.0 * prop(:,nj-1) & 
                           + 1022.0 * prop(:,nj-2) - 258.0 * prop(:,nj-3) &
                           + 77.0 * prop(:,nj-4) - 11.0 * prop(:,nj-5) ) / 1440.0
          flux_i(:,1) = ( 475.0 * prop(:,1) + 1427.0 * prop(:,2) &
                        - 798.0 * prop(:,3) + 482.0 * prop(:,4) &
                        - 173.0 * prop(:,5) + 27.0 * prop(:,6) ) / 1440.0
          flux_i(:,nj-1) = ( 475.0 * prop(:,nj) + 1427.0 * prop(:,nj-1) &
                           - 798.0 * prop(:,nj-2) + 482.0 * prop(:,nj-3) &
                           - 173.0 * prop(:,nj-4) + 27.0 * prop(:,nj-5) ) / 1440.0     
                    
      else if(av%spatial_acc == 4) then
          
          flux_i(:,2:nj-2) = ( - prop(:,1:nj-3) + 13.0 * prop(:,2:nj-2) &
                             + 13.0 * prop(:,3:nj-1) - prop(:,4:nj)) / 24.0
          flux_i(:,1) = ( 9.0 * prop(:,1) + 19.0 * prop(:,2) &
                        - 5.0 * prop(:,3) + prop(:,4) ) / 24.0
          flux_i(:,nj-1) = ( 9.0 * prop(:,nj) + 19.0 * prop(:,nj-1) &
                           - 5.0 * prop(:,nj-2) + prop(:,nj-3) ) / 24.0  
          
      else !if no valid spatial accuracy set then use 2nd order
          flux_i(:,1:nj-1) = ( prop(:,1:nj-1) + prop(:,2:nj) ) / 2.0
      end if
      
      end subroutine get_flux_i
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      subroutine get_flux_j(av,prop,flux_j)
      
!     gets the i flux component using either fourth or second order estimate

      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(in) :: prop(:,:)
      real, intent(inout) :: flux_j(:,:)
      integer :: ni, nj
      
      ni = size(prop,1); nj = size(prop,2)
      
      if(av%spatial_acc == 6) then
          
          flux_j(3:ni-3,:) = ( 11.0 * prop(1:ni-5,:) - 93.0 * prop(2:ni-4,:) & 
                             + 802.0 * ( prop(3:ni-3,:) + prop(4:ni-2,:) ) &
                             - 93.0 * prop(5:ni-1,:) + 11.0 * prop(6:ni,:) ) / 1440.0
          flux_j(2,:) = ( - 27.0 * prop(1,:) + 637.0 * prop(2,:) & 
                        + 1022.0 * prop(3,:) - 258.0 * prop(4,:) &
                        + 77.0 * prop(5,:) - 11.0 * prop(6,:) ) / 1440.0
          flux_j(ni-2,:) = ( - 27.0 * prop(ni,:) + 637.0 * prop(ni-1,:) & 
                           + 1022.0 * prop(ni-2,:) - 258.0 * prop(ni-3,:) &
                           + 77.0 * prop(ni-4,:) - 11.0 * prop(ni-5,:) ) / 1440.0
          flux_j(1,:) = ( 475.0 * prop(1,:) + 1427.0 * prop(2,:) &
                        - 798.0 * prop(3,:) + 482.0 * prop(4,:) &
                        - 173.0 * prop(5,:) + 27.0 * prop(6,:) ) / 1440.0
          flux_j(ni-1,:) = ( 475.0 * prop(ni,:) + 1427.0 * prop(ni-1,:) &
                           - 798.0 * prop(ni-2,:) + 482.0 * prop(ni-3,:) &
                           - 173.0 * prop(ni-4,:) + 27.0 * prop(ni-5,:) ) / 1440.0            
                
      else if(av%spatial_acc == 4) then
          
          flux_j(2:ni-2,:) = (- prop(1:ni-3,:) + 13.0 * prop(2:ni-2,:) &
                             + 13.0 * prop(3:ni-1,:) - prop(4:ni,:)) / 24.0
          flux_j(1,:) = (9.0 * prop(1,:) + 19.0 * prop(2,:) &
                        - 5.0 * prop(3,:) + prop(4,:) ) / 24.0
          flux_j(ni-1,:) = (9.0 * prop(ni,:) + 19.0 * prop(ni-1,:) &
                           - 5.0 * prop(ni-2,:) + prop(ni-3,:) ) / 24.0 
           
      else !if no valid spatial accuracy set then use 2nd order
          flux_j(1:ni-1,:) = ( prop(1:ni-1,:) + prop(2:ni,:) ) / 2.0
      end if
      
      end subroutine get_flux_j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module flux_stencil


