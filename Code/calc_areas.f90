      
      subroutine calc_areas(g)

!     Calculate the area of the quadrilateral cells and the lengths of the side
!     facets

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_grid), intent(inout) :: g
      integer :: ni, nj

!     Declare integers or any extra variables you need here
      

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Calculate the areas of the cells and store in g%area. The area of any
!     quadrilateral is half of the magnitude of the cross product of the two
!     vectors that form the diagonals. Check the order of your product so that
!     the values come out positive! You can do this using two nested loops in
!     the i and j-directions or in a vectorised way by indexing the coordinate
!     arrays with lists of indices
      g%area = 0.5 * ( (g%x(2:ni,2:nj) - g%x(1:ni-1,1:nj-1)) &
                         * (g%y(1:ni-1,2:nj) - g%y(2:ni,1:nj-1)) &
                        - ((g%x(1:ni-1,2:nj) - g%x(2:ni,1:nj-1)) &
                         * (g%y(2:ni,2:nj) - g%y(1:ni-1,1:nj-1))))

!     Calculate the projected lengths in the x and y-directions on all of the
!     "i = const" facets and store them in g%lx_i and g%ly_i. When combined
!     together these two components define a vector that is normal to the facet,
!     pointing inwards towards the centre of the cell. This is only the case for
!     the left hand side of the cell, the vector stored in position i,j points
!     towards the centre of the i,j cell
      g%lx_i = g%x(1:ni-1,2:nj) - g%x(1:ni-1,1:nj-1)
      g%ly_i = g%y(1:ni-1,2:nj) - g%y(1:ni-1,1:nj-1)

!     Now repeat the calculation for the project lengths on the "j=const"
!     facets. 
      g%lx_j = g%x(2:ni,1:nj-1) - g%x(1:ni-1,1:nj-1)
      g%ly_j = g%y(2:ni,1:nj-1) - g%y(1:ni-1,1:nj-1)

!     Find the minimum length scale in the mesh, this is defined as the length
!     of the shortest side of all the cells. Call this length "l_min", it is used
!     to set the timestep from the CFL number. Start by calculating the lengths
!     of the i and j facets by using the intrinsic function "hypot", this avoids
!     underflow and overflow errors. Then find the overal minimum value using
!     both the "min" and "minval" functions.
      l_i = hypot(g%lx_i(:,:),g%ly_i(:,:))
      l_j = hypot(g%lx_j(:,:),g%ly_j(:,:))
      
      g%l_min = min(minval(l_i),minval(l_j))
!
!     Print the overall minimum length size that has been calculated
      write(6,*) 'Calculated cell areas and facet lengths'
      write(6,*) '  Overall minimum element size = ', g%l_min
      write(6,*)

      end subroutine calc_areas
