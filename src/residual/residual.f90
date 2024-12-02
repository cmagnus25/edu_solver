!********************************************************************************

!********************************************************************************
!* This subroutine computes the residual for a node-centered finite-volume method
!*
!* ------------------------------------------------------------------------------
!*  Input: the current solution
!*
!* Output: node(:)%res = the residual computed by the current solution.
!* ------------------------------------------------------------------------------
!*
!* Note: dU/dt + dF/dx + dG/dy = 0. Residuals are first computed as
!*       the integral of (dF/dx + dG/dy), and at the end negative sign is added
!*       so that we have dU/dt = Res at every node.
!* 
!********************************************************************************
module residual

 private
 
 public :: compute_residual_ncfv
 
 contains

 subroutine compute_residual_ncfv

 use edu2d_constants   , only : p2, zero, half
 use edu2d_my_main_data, only : nnodes, node, nedges, edge, nbound, bound, &
                                rho_inf, u_inf, v_inf, p_inf, nq 
 use input_parameter, only : inviscid_flux, gradient_type 						  
 use gradients_lsq  , only : compute_gradient_nc
 use roe_scheme     , only : roe
 
 implicit none

!Local variables
 real(p2), dimension(nq) :: num_flux       !Numerical flux
 real(p2), dimension(nq) :: wL, wR         !Left and right face values
 real(p2), dimension(nq) :: dwL, dwR       !Slope at left and right nodes
 real(p2), dimension(2)  :: e12            !Unit edge vector
 real(p2), dimension(2)  :: n12            !Unit face normal vector
 real(p2)                :: mag_e12        !Magnitude of the edge vector
 real(p2)                :: mag_n12        !Magnitude of the face-normal vector
 real(p2)                :: wsn            !Scaled maximum wave speed
 real(p2), dimension(nq) :: bfluxL, bfluxR !Boundary flux at left/right nodes
 integer                 :: node1, node2   !Left and right nodes of each edge
 integer                 :: n1, n2         !Left and right nodes of boundary face
 integer                 :: i, j
 integer                 :: ix=1, iy=2

!-------------------------------------------------------------------------
! Gradient Reconstruction for second-order accuracy

!  Initialization

    nodes0 : do i = 1, nnodes
     node(i)%gradw = zero
    end do nodes0

!  Perform LSQ gradient computations in the premitive variables w=[rho,u,v,p]:

     call compute_gradient_nc(1,gradient_type) ! Density gradient: grad(rho)
     call compute_gradient_nc(2,gradient_type) !Velocity gradient: grad(u)
     call compute_gradient_nc(3,gradient_type) !Velocity gradient: grad(v)
     call compute_gradient_nc(4,gradient_type) !Pressure gradient: grad(p)

!-------------------------------------------------------------------------
! Residual computation: interior fluxes

! Initialization
  nodes : do i = 1, nnodes
   node(i)%res = zero
   node(i)%wsn = zero
  end do nodes

! Flux computation across internal edges (to be accumulated in res(:))
!
!   node2              1. Extrapolate the solutions to the edge-midpoint
!       o                 from the nodes, n1 and n2.
!        \   face      2. Compute the numerical flux
!         \ -------c2  3. Add it to the residual for n1, and subtract it from
!        / \              the residual for n2.
!   face/   \ edge
!      /     o         Directed area is the sum of the left and the right faces.
!    c1    node1       Left/right face is defined by the edge-midpoint and
!                      the centroid of the left/right element.
!                      Directed area is positive in n1 -> n2
!
! (c1, c2: element centroids)
!
!--------------------------------------------------------------------------------
  edges : do i = 1, nedges
!--------------------------------------------------------------------------------

! Left and right nodes of the i-th edge

    node1 = edge(i)%n1  ! Left node of the edge
    node2 = edge(i)%n2  ! Right node of the edge
      n12 = edge(i)%dav ! This is the directed area vector (unit vector)
  mag_n12 = edge(i)%da  ! Magnitude of the directed area vector
      e12 = edge(i)%ev  ! This is the vector along the edge (uniti vector)
  mag_e12 = edge(i)%e   ! Magnitude of the edge vector (Length of the edge)

! Solution gradient projected along the edge
!
!  NOTE: The gradient is multiplied by the distance.
!        So, it is equivalent to the solution difference.

  dwL = (node(node1)%gradw(:,ix)*e12(ix) + node(node1)%gradw(:,iy)*e12(iy) )*half*mag_e12
  dwR = (node(node2)%gradw(:,ix)*e12(ix) + node(node2)%gradw(:,iy)*e12(iy) )*half*mag_e12

!  No limiter (good for smooth solutions); You can implement a limiter here.

       wL = node(node1)%w + dwL
       wR = node(node2)%w - dwR

!  Compute the numerical flux for given wL and wR.

!  (1) Roe flux (carbuncle is expected for strong shocks)
   if     (trim(inviscid_flux)=="roe") then

     call roe(wL,wR,n12, num_flux,wsn)

   else

    write(*,*) " Invalid input for inviscid_flux = ", trim(inviscid_flux)
    write(*,*) " Choose roe (no others available), and try again."
    write(*,*) " ... Stop."
    stop

   endif

!  Add the flux multiplied by the magnitude of the directed area vector to node1,
!  and accumulate the max wave speed quantity 'wsn' for use in the time step calculation.

     node(node1)%res = node(node1)%res  +  num_flux * mag_n12
     node(node1)%wsn = node(node1)%wsn  +       wsn * mag_n12

!  Subtract the flux multiplied by the magnitude of the directed area vector from node2,
!  and accumulate the max wave speed quantity 'wsn' for use in the time step calculation.
!
!  NOTE: Subtract because the outward face normal is -n12 for the node2.

     node(node2)%res = node(node2)%res  -  num_flux * mag_n12
     node(node2)%wsn = node(node2)%wsn  +       wsn * mag_n12

!--------------------------------------------------------------------------------
  end do edges
!--------------------------------------------------------------------------------

!-------------------------------------------------------------------------
! Close with the boundary flux using the element-based formula that is
! exact for linear fluxes (See Nishikawa AIAA2010-5093 for boundary weights
! that ensure the linear exactness for 2D/3D elements).
! Here, the formula for triangles is implemented.
!
!
!      |  Interior Domain          |
!      |        .........          |
!      |        .       .          |
!      |        .       .          |
!      o--o--o-----o---------o--o--o  <- Boundary segment
!                  n1   |   n2
!                       v
!                     n12 (unit face normal vector)
!
! NOTE: We visit each boundary face, defined by the nodes n1 and n2,
!       and compute the flux across the boundary face: left half for node1,
!       and the right half for node2. In the above figure, the dots indicate
!       the control volume around the node n1. Clearly, the flux across the
!       left half of the face contributes to the node n1. Similarly for n2.
!
!
!--------------------------------------------------------------------------------
  bc_loop : do i = 1, nbound
!--------------------------------------------------------------------------------

!------------------------------------------------
!  BC: Upwind flux via free stream values
!
!      NOTE: If the final solution at the boundary node is far from
!            the free stream values, then the domain is probably is not large enough.

   bc : if (trim(bound(i)%bc_type) == "freestream") then

    bnodes_numerical_flux_via_freestream : do j = 1, bound(i)%nbfaces

         n1 = bound(i)%bnode(j  )  !Left node
         n2 = bound(i)%bnode(j+1)  !Right node
     n12(1) = bound(i)%bfnx(j)     !x-component of the unit face normal vector
     n12(2) = bound(i)%bfny(j)     !y-component of the unit face normal vector
    mag_n12 = bound(i)%bfn(j)*half !Half length of the boundary face, j.

!   1. Left node
       wL = node(n1)%w
       wR = (/ rho_inf, u_inf, v_inf, p_inf /) ! = w_out, state outside
       call roe(wL,wR,n12, num_flux,wsn)
       bfluxL = num_flux
       node(n1)%wsn = node(n1)%wsn + wsn*mag_n12

!   2. Right node
       wL = node(n2)%w
       wR = (/ rho_inf, u_inf, v_inf, p_inf /) ! = w_out, state outside
       call roe(wL,wR,n12, num_flux,wsn)
       bfluxR = num_flux
       node(n2)%wsn = node(n2)%wsn + wsn*mag_n12

!   3. Add contributions to the two nodes.
!      Note: This is the formula for triangles, which is exact for linear fluxes.
!            For other elements, suitable formulas can be found in Nishikawa AIAA2010-5093.

       node(n1)%res = node(n1)%res + (5.0_p2*bfluxL + bfluxR)/6.0_p2*mag_n12
       node(n2)%res = node(n2)%res + (5.0_p2*bfluxR + bfluxL)/6.0_p2*mag_n12

    end do bnodes_numerical_flux_via_freestream

!------------------------------------------------
!  BC: Weak condition on a slip solid wall
!
!      NOTE: Physical flux has only pressure flux since qn=0.

   elseif (trim(bound(i)%bc_type) == "slip_wall_weak") then

    bnodes_slip_wall_weak : do j = 1, bound(i)%nbfaces

         n1 = bound(i)%bnode(j  ) ! Left node
         n2 = bound(i)%bnode(j+1) ! Right node
     n12(1) = bound(i)%bfnx(j)
     n12(2) = bound(i)%bfny(j)
    mag_n12 = bound(i)%bfn(j)*half

!   1. Left node

       num_flux(1)  = zero
       num_flux(2)  = node(n1)%w(4)*n12(1) !Pressure flux only
       num_flux(3)  = node(n1)%w(4)*n12(2) !Pressure flux only
       num_flux(4)  = zero

             bfluxL = num_flux

                wsn = sqrt(1.4_p2*node(n1)%w(4)/node(n1)%w(1)) !Speed of sound
       node(n1)%wsn = node(n1)%wsn + wsn*mag_n12

!   2. Right node

       num_flux(1)  = zero
       num_flux(2)  = node(n2)%w(4)*n12(1) !Pressure flux only
       num_flux(3)  = node(n2)%w(4)*n12(2) !Pressure flux only
       num_flux(4)  = zero

             bfluxR = num_flux

                wsn = sqrt(1.4_p2*node(n2)%w(4)/node(n2)%w(1)) !Speed of sound
       node(n2)%wsn = node(n2)%wsn + wsn*mag_n12

!   3. Add contributions to the two nodes.
!      Note: This is the formula for triangles, which is exact for linear fluxes.
!            For other elements, suitable formulas can be found in Nishikawa AIAA2010-5093.

       node(n1)%res = node(n1)%res + (5.0_p2*bfluxL + bfluxR)/6.0_p2*mag_n12
       node(n2)%res = node(n2)%res + (5.0_p2*bfluxR + bfluxL)/6.0_p2*mag_n12

    end do bnodes_slip_wall_weak

!------------------------------------------------

   endif bc

!--------------------------------------------------------------------------------
  end do bc_loop
!--------------------------------------------------------------------------------

! Switch the residual sign, which goes to the right hand side of a linear system.

  nodes3 : do i = 1, nnodes

   node(i)%res = - node(i)%res

  end do nodes3

 end subroutine compute_residual_ncfv
!--------------------------------------------------------------------------------
end module residual