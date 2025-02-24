
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
      real, intent(out) :: dcell(:,:)
      real, dimension(size(prop,1),size(prop,2)) :: dnode
      integer :: ni, nj

!     Get the block size and store locally for convenience
      ni = size(prop,1); nj = size(prop,2)

!     Use the finite volume method to find the change in the variables "prop"
!     over the timestep "dt", save it in the array "dcell"
      dcell = ( av%dt / area ) * & 
         (flux_i(1:ni-1,:) - flux_i(2:ni,:) + flux_j(:,1:nj-1) - flux_j(:,2:nj))

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

      end subroutine sum_fluxes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module flux_stencil


