!********************************************************************************
!* Educationally-Designed Unstructured 2D (EDU2D) Code
!*
!*
!*     This module belongs to the inviscid version: EDU2D-Euler-Implicit
!*
!*
!*
!* This module contains all subroutines needed to advance the solution in time.
!*
!* Contained subroutines/functions:
!* -----------------------------------------------------------------------------
!* euler_solver_main          : Main time-stepping routine
!* compute_residual_ncfv      : Compute the residual (RHS, spatial part)
!* update_solution            : Update the solution at all nodes
!* update_solution_du         : Update the solution at all nodes by precomputed du
!* roe                        : Roe flux with an entropy fix
!* compute_time_step          : Compute global/local time step
!* residual_norm              : Compute the residual norms (L1, L2, Linf)
!* w2u                        : Primitive to conservative variables
!* u2w                        : Conservative to primitive variables
!*
!* And subroutines for LSQ gradiets
!* -----------------------------------------------------------------------------
!*
!* Written by Katate Masatsuka (http://www.cfdbooks.com)
!********************************************************************************
 module edu2d_euler_implct_solver

 private

 public :: euler_solver_main
 public :: compute_lsq_coeff_nc
 public ::   check_lsq_coeff_nc
 public :: w2u

 contains

!********************************************************************************
!********************************************************************************
!********************************************************************************
!* Euler solver: Node-Centered Finite-Volume Method (Edge-Based)
!*
!* - Node-centered second-order finite-volume method for unstructured grids
!* - Roe flux with an entropy fix
!* - Reconstruction by unweighted least-squares method (2x2 system for gradients)
!* - Explicit: 2-Stage Runge-Kutta time-stepping
!* - Implicit: Defect correction with a pseudo time term
!*
!********************************************************************************
 subroutine euler_solver_main

 use edu2d_constants   , only : p2, one, half

 use edu2d_my_main_data, only : CFL, nnodes, node, tolerance, cfl1, cfl2, CFL_ramp_steps, &
                                iteration_method, max_iterations, CFLexp, nq

!Implicit method uses subroutines below to construct Jacobian matrix
!and relax the linear system.
 use edu2d_euler_jacobian    , only : construct_jacobian_ncfv ! Construct Jacobian
 use edu2d_euler_linear_solve, only : gs_sequential           ! GS relaxaton for linear system

 implicit none

!Local variables
 real(p2), dimension(nq,3)             :: res_norm      !Residual norms(L1,L2,Linf)
 real(p2), dimension(:,:), allocatable :: u0            !Saved solution
 integer                               :: i_iteration   !Iteration counter
 real(p2)                              :: s, exp_factor !CFL ramping variables
 integer                               :: i
 real(p2)                              :: time_begin    !Starting time
 real(p2)                              :: time_end      !End time

! For explicit scheme:
!  1. Allocate the temporary solution array needed for the Runge-Kutta method.
!  2. Assign the CFL number = CFLexp (which is fixed and never changes).

  if (trim(iteration_method) == "explicit") then

    allocate(u0(nnodes,nq))
    CFL = CFLexp

  endif

!--------------------------------------------------------------------------------
! Iteration towards steady state
!--------------------------------------------------------------------------------
  call cpu_time( time_begin )

  iteration : do i_iteration = 1, max_iterations

!*********************************************************************
!*********************************************************************
! Method 1: Explicit time stepping (RK2) with local time stepping

  method : if (trim(iteration_method) == "explicit") then

  !------------------------------------------------------
  ! Two-stage Runge-Kutta scheme: u^n is saved as u0(:,:)
  !  1. u^*     = u^n - (dt/vol)*Res(u^n)
  !  2. u^{n+1} = 1/2*u^n + 1/2*[u^* - (dt/vol)*Res(u^*)]
  !------------------------------------------------------

  !-----------------------------
  !- 1st Stage of Runge-Kutta:
  !-----------------------------

!   Compute Res(u^n)
    call compute_residual_ncfv

!   Compute residual norms (undivided residual)
    call residual_norm(res_norm)

!   Display the residual norm.
    if (i_iteration==1) then
      write(*,'(a94)') "Density    X-momentum  Y-momentum   Energy"
    endif

    if (mod(i_iteration-1,1)==0) write(*,'(a4,es12.2,a13,i6,a12,4es12.2)') &
       "CFL=", CFLexp, "Iteration=", i_iteration-1, " L1(res)=",res_norm(:,1)

!   Stop if the L1 residual norm drops below the specified tolerance
    if (maxval(res_norm(:,1)) < tolerance) exit iteration

!   Save off the previous solution (to be used in the 2nd stage).
    do i = 1, nnodes
     u0(i,:) = node(i)%u
    end do

!   Compute the time step (local=node(i)%dt)
    call compute_time_step

!   Update the solution by the local time step
!   1st Stage => u^* = u^n - dt/dx*Res(u^n)
    call update_solution(one)

  !-----------------------------
  !- 2nd Stage of Runge-Kutta:
  !-----------------------------

!   Compute Res(u^*)
    call compute_residual_ncfv

!   Compute 1/2*(u^n + u^*)
    do i = 1, nnodes
     node(i)%u = half*( node(i)%u + u0(i,:) )
    end do

!   2nd Stage => u^{n+1} = 1/2*(u^n + u^*) - 1/2*dt/dx*Res(u^*)
    call update_solution(half)

!   Go to the next iteration (local time step)
    cycle iteration

!*********************************************************************
!*********************************************************************
! Method 2: Implicit (Defect Correction - DC, with a pseudo time term)
!
  elseif (trim(iteration_method) == "implicit") then

!   Exponential CFL ramping between CFL1 and CFL2 over 'CFL_ramp_steps' iterations.
!   Note: This is by no means an optimal ramping strategy.
    if (i_iteration <= CFL_ramp_steps) then

     exp_factor = 0.5_p2
       s = real( i_iteration-1, p2) / real( CFL_ramp_steps-1, p2)
     CFL = CFL1 + (CFL2-CFL1)*( one - exp(-s*exp_factor) )/ ( one - exp(-exp_factor) )

    endif

!   Compute Residual: Right hand side of [V/dt+dR/dU]*dU = -Residual
    call compute_residual_ncfv   ! -> This computes node(:)%res.

!   Compute residual norms (undivided residual)
    call residual_norm(res_norm) ! -> This computes res_norm.

!   Display the residual norm.
    if (i_iteration==1) then
      write(*,'(a101)') "Density    X-momentum  Y-momentum   Energy"
    endif

    if (mod(i_iteration-1,1)==0) write(*,'(a4,f20.2,a13,i6,a12,4es12.2)') &
       "CFL=", CFL, "iteration=", i_iteration-1, " L1(res)=",res_norm(:,1)

!   Stop if the L1 residual norm drops below the specified tolerance for all eqns.
    if (maxval(res_norm(:,1)) < tolerance) exit iteration

!   Compute the local time step (the local time step is used by DC in pseudo time term)
    call compute_time_step       ! -> This computes node(:)%dt.

!   Construct the residual Jacobian matrix (based on 1st-order scheme: Defect Correction)
!   which is the matrix on the left hand side of [V/dt+dR/dU]*dU = -Residual.
!   Note: This sbroutine can be found in edu2d_euler_jacobian_v0.f90
    call construct_jacobian_ncfv ! -> This computes jac(:)%diag and jac(:)%off.

!   Relax the linear system by sequential Gauss-Seidel to get du (Correction)
!   Note: This subroutine can be found in edu2d_euler_inear_solve_v0.f90
    call gs_sequential           ! -> This relaxes the linear system and compute du.

!   Update the solution: u_new = u_old + du
    call update_solution_du      ! -> This updates node(:)%u.

!   Go to the next iteration
    cycle iteration

  endif method

  end do iteration

!--------------------------------------------------------------------------------
! End of iteration toward steady state
!--------------------------------------------------------------------------------

  write(*,*)
  call cpu_time( time_end )
  write(*,*) " Total CPU time to solution = ", time_end - time_begin, " seconds"
  write(*,*)

  if (i_iteration == max_iterations) then
   write(*,*) " Not converged... Sorry..."
   write(*,*) "   Max iterations reached... max_iterations=", max_iterations
   write(*,*) "   Increase max_iterations, and try again."
  endif

  write(*,*) " Converged."
  write(*,*) " Final iteration      =", i_iteration-1
  write(*,*)
  write(*,*) "Finished the Euler solver... Bye!"

 end subroutine euler_solver_main
!--------------------------------------------------------------------------------
!********************************************************************************
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
 subroutine compute_residual_ncfv

 use edu2d_constants   , only : p2, zero, half
 use edu2d_my_main_data, only : nnodes, node, nedges, edge, nbound, bound, &
                          rho_inf, u_inf, v_inf, p_inf, &
                          inviscid_flux, nq, gradient_type
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


!********************************************************************************
!* This subroutine updates the solution for explicit scheme (RK2).
!*
!*  Note: This is not used by implicit scheme.
!*
!* ------------------------------------------------------------------------------
!*  Input:       coeff = coefficient for RK time-stepping
!*                 CFL = CFL number
!*          node(:)res = the residual
!*
!* Output:   node(:)u  = updated conservative variables 
!*           node(:)w  = updated primitive variables 
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine update_solution(coeff)

 use edu2d_constants   , only : p2, zero
 use edu2d_my_main_data, only : nnodes, node, M_inf

 implicit none

 real(kind=p2), intent(in) :: coeff

!Local variables
 integer :: i

!--------------------------------------------------------------------------------
  nodes : do i = 1, nnodes
!--------------------------------------------------------------------------------

!   Solution change based on the local time step, node(i)%dt.

     node(i)%du = (coeff*node(i)%dt/node(i)%vol) * node(i)%res

!   Note: The explicit scheme is not time-accurate because the time step
!         is not global. For time-accurate computation, use the following:
!    node(i)%du = (coeff*dt/node(i)%vol) * node(i)%res

!   Solution update

     node(i)%u  = node(i)%u + node(i)%du !Make changes
     node(i)%w  = u2w(node(i)%u)         !Update primitive variables

!--------------------------------------------------------------------------------
  end do nodes
!--------------------------------------------------------------------------------

! Check the density/pressure positivity.
!

  nodes2 : do i = 1, nnodes

!  Check density
   if (node(i)%w(1) <= zero) then
    write(*,*) " Negative density detected in update(): rho = ", &
               node(i)%w(1), i, nnodes
    node(i)%w(1) = min(0.1_p2*M_inf**2, 1.0e-04_p2)
      node(i)%u  = w2u(node(i)%w)
    write(*,*) " Stopped at 'update_solution()'"
    stop
   endif

!  Check pressure
   if (node(i)%w(4) <= zero) then
    write(*,*) " Negative pressure detected in update(): p = ", &
               node(i)%w(4), i, nnodes
    node(i)%w(4) = min(0.1_p2*M_inf**2, 1.0e-04_p2)
      node(i)%u  = w2u(node(i)%w)
    write(*,*) " Stopped at 'update_solution()'"
    stop
   endif

  end do nodes2


 end subroutine update_solution
!--------------------------------------------------------------------------------


!********************************************************************************
!* This subroutine updates the solution for implicit scheme (DC).
!*
!*  Note: This is not used by explicit scheme.
!*
!* ------------------------------------------------------------------------------
!*  Input:  node(:)du  = correction computed by relaxing the linear system
!*
!* Output:   node(:)u  = updated conservative variables 
!*           node(:)w  = updated primitive variables 
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine update_solution_du

 use edu2d_constants   , only : p2, zero
 use edu2d_my_main_data, only : nnodes, node, M_inf

 implicit none

!Local variables
 integer                :: i

!--------------------------------------------------------------------------------
  nodes : do i = 1, nnodes
!--------------------------------------------------------------------------------

!   Solution update

     node(i)%u  = node(i)%u + node(i)%du !Make changes by du computed by DC
     node(i)%w  = u2w(node(i)%u)         !Update primitive variables

!--------------------------------------------------------------------------------
  end do nodes
!--------------------------------------------------------------------------------

! Check the density/pressure positivity.
!

  nodes2 : do i = 1, nnodes

!  Check density
   if (node(i)%w(1) <= zero) then
    write(*,*) " Negative density detected in update_solution_du(): rho = ", &
               node(i)%w(1), i, nnodes
    node(i)%w(1) = min(0.1_p2*M_inf**2, 1.0e-04_p2)
      node(i)%u  = w2u(node(i)%w)
    write(*,*) " Stopped at 'update_solution_du()'"
    stop
   endif

!  Check pressure
   if (node(i)%w(4) <= zero) then
    write(*,*) " Negative pressure detected in update_solution_du(): p = ", &
               node(i)%w(4), i, nnodes
    node(i)%w(4) = min(0.1_p2*M_inf**2, 1.0e-04_p2)
      node(i)%u  = w2u(node(i)%w)
    write(*,*) " Stopped at 'update_solution_du()'"
    stop
   endif

  end do nodes2


 end subroutine update_solution_du
!--------------------------------------------------------------------------------


!********************************************************************************
!* -- Roe's Flux Function with entropy fix---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!*
!* NOTE: 3D version of this subroutine is available for download at
!*       http://cfdbooks.com/cfdcodes.html
!*
!* ------------------------------------------------------------------------------
!*  Input:   primL(1:4) =  left state (rhoL, uL, vL, pL)
!*           primR(1:4) = right state (rhoR, uR, vR, pR)
!*               njk(2) = Face normal (L -> R). Must be a unit vector.
!*
!* Output:    flux(1:4) = numerical flux
!*                  wsn = half the max wave speed
!*                        (to be used for time step calculations)
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine roe(primL, primR, njk,  flux, wsn)

 use edu2d_constants   , only : p2, zero, one, two, half, fifth
 use edu2d_my_main_data, only : gamma

 implicit none

!Input:
 real(p2), intent( in) :: primL(4), primR(4) ! [rho, u, v, p]_{L,R}
 real(p2), intent( in) :: njk(2)             ! Face normal, njk=[nx, ny]

!Output
 real(p2), intent(out) :: flux(4)
 real(p2), intent(out) :: wsn

!Local variables
 real(p2) :: nx, ny                  ! Normal vector
 real(p2) :: mx, my                  ! Tangent vector: mx*nx+my*ny = 0
 real(p2) :: uL, uR, vL, vR          ! Velocity components.
 real(p2) :: rhoL, rhoR, pL, pR      ! Primitive variables.
 real(p2) :: unL, unR, umL, umR      ! Normal and tangent velocities
 real(p2) :: aL, aR, HL, HR          ! Speeds of sound.
 real(p2) :: RT,rho,u,v,H,a,un, um   ! Roe-averages
 real(p2) :: drho,dun,dum,dp,LdU(4)  ! Wave strenghs
 real(p2) :: ws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
 real(p2) :: fL(4), fR(4), diss(4)   ! Fluxes ad dissipation term
 real(p2) :: dws(4)                  ! User-specified width for entropy fix
 integer :: i, j

  nx = njk(1)
  ny = njk(2)

!Tangent vector (Do you like it? Actually, Roe flux can be implemented 
! without any tangent vector. See "I do like CFD, VOL.1" for details.)
  mx = -ny
  my =  nx

!Primitive and other variables.
!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
     unL = uL*nx+vL*ny
     umL = uL*mx+vL*my
      pL = primL(4)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL)
!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
     unR = uR*nx+vR*ny
     umR = uR*mx+vR*my
      pR = primR(4)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR)

!First compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
     u = (uL+RT*uR)/(one+RT)
     v = (vL+RT*vR)/(one+RT)
     H = (HL+RT* HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*(u*u+v*v)) )
    un = u*nx+v*ny
    um = u*mx+v*my

!Wave Strengths
   drho = rhoR - rhoL 
     dp =   pR - pL
    dun =  unR - unL
    dum =  umR - umL

  LdU(1) = (dp - rho*a*dun )/(two*a*a)
  LdU(2) = rho*dum
  LdU(3) =  drho - dp/(a*a)
  LdU(4) = (dp + rho*a*dun )/(two*a*a)

!Wave Speed
  ws(1) = abs(un-a)
  ws(2) = abs(un)
  ws(3) = abs(un)
  ws(4) = abs(un+a)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
! only for the nonlinear fields.
  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(4) = fifth
   if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

!Right Eigenvectors
  Rv(1,1) = one    
  Rv(2,1) = u - a*nx
  Rv(3,1) = v - a*ny
  Rv(4,1) = H - un*a

  Rv(1,2) = zero
  Rv(2,2) = mx
  Rv(3,2) = my
  Rv(4,2) = um

  Rv(1,3) = one
  Rv(2,3) = u
  Rv(3,3) = v 
  Rv(4,3) = half*(u*u+v*v)

  Rv(1,4) = one
  Rv(2,4) = u + a*nx
  Rv(3,4) = v + a*ny
  Rv(4,4) = H + un*a

!Dissipation Term
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*LdU(j)*Rv(i,j)
   end do
  end do

!Compute the flux.
  fL(1) = rhoL*unL
  fL(2) = rhoL*unL * uL + pL*nx
  fL(3) = rhoL*unL * vL + pL*ny
  fL(4) = rhoL*unL * HL

  fR(1) = rhoR*unR
  fR(2) = rhoR*unR * uR + pR*nx
  fR(3) = rhoR*unR * vR + pR*ny
  fR(4) = rhoR*unR * HR

  flux = half * (fL + fR - diss)
  wsn = half*(abs(un) + a)  !Normal max wave speed times half

 end subroutine roe
!--------------------------------------------------------------------------------

!********************************************************************************
!* This subroutine computes the local time-step at all nodes.
!*
!* ------------------------------------------------------------------------------
!*  Input: node(i)%vol = Dual volume
!*         node(i)%wsn = Sum of the max wave speed multiplied by the face length
!*
!* Output: node(:)%dt  =  local time step
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine compute_time_step

 use edu2d_my_main_data, only : nnodes, node, CFL

 implicit none

!Local variables
 integer  :: i

!--------------------------------------------------------------------------------
  nodes : do i = 1, nnodes
!--------------------------------------------------------------------------------

! Local time step: dt at node i = volume/sum(0.5*max_wave_speed*face_area).

    node(i)%dt = CFL*node(i)%vol / node(i)%wsn

!--------------------------------------------------------------------------------
  end do nodes
!--------------------------------------------------------------------------------

 end subroutine compute_time_step
!--------------------------------------------------------------------------------

!********************************************************************************
!* This subroutine computes the residual norms: L1, L2, L_infty
!*
!* ------------------------------------------------------------------------------
!*  Input:  node(:)res = the residuals
!*
!* Output:  res_norm   = residual norms (L1, L2, Linf)
!* ------------------------------------------------------------------------------
!*
!* NOTE: It is not done here, but I advise you to keep the location of the
!*       maximum residual (L_inf).
!*
!********************************************************************************
 subroutine residual_norm(res_norm)

 use edu2d_constants   , only : p2, zero, one
 use edu2d_my_main_data, only : node, nnodes

 implicit none

 real(kind=p2), dimension(4,3), intent(out)   :: res_norm

!Local variables
 real(kind=p2), dimension(4) :: residual
 integer :: i

   res_norm(:,1) =  zero
   res_norm(:,2) =  zero
   res_norm(:,3) = - one

!--------------------------------------------------------------------------------
  nodes : do i = 1, nnodes
!--------------------------------------------------------------------------------
   residual = abs( node(i)%res )                  !Univided residual
   res_norm(:,1) = res_norm(:,1)    + residual    !L1   norm
   res_norm(:,2) = res_norm(:,2)    + residual**2 !L2   norm
   res_norm(:,3) = max(res_norm(:,3), residual)   !Linf norm
!--------------------------------------------------------------------------------
  end do nodes
!--------------------------------------------------------------------------------

   res_norm(:,1) =      res_norm(:,1)/real(nnodes,p2)
   res_norm(:,2) = sqrt(res_norm(:,2)/real(nnodes,p2))

 end subroutine residual_norm
!--------------------------------------------------------------------------------


!********************************************************************************
!* Compute U from W
!*
!* ------------------------------------------------------------------------------
!*  Input:  w =    primitive variables (rho,     u,     v,     p)
!* Output:  u = conservative variables (rho, rho*u, rho*v, rho*E)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
 function w2u(w) result(u)

 use edu2d_constants   , only : p2, one, half
 use edu2d_my_main_data, only : gamma

 implicit none

 real(p2), dimension(4), intent(in) :: w

!Local variables
 real(p2), dimension(4)             :: u !output

  u(1) = w(1)
  u(2) = w(1)*w(2)
  u(3) = w(1)*w(3)
  u(4) = w(4)/(gamma-one)+half*w(1)*(w(2)*w(2)+w(3)*w(3))

 end function w2u
!--------------------------------------------------------------------------------

!********************************************************************************
!* Compute U from W
!*
!* ------------------------------------------------------------------------------
!*  Input:  u = conservative variables (rho, rho*u, rho*v, rho*E)
!* Output:  w =    primitive variables (rho,     u,     v,     p)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
 function u2w(u) result(w)

 use edu2d_constants   , only : p2, one, half
 use edu2d_my_main_data, only : gamma

 implicit none

 real(p2), dimension(4), intent(in) :: u

!Local variables
 real(p2), dimension(4)             :: w !output

  w(1) = u(1)
  w(2) = u(2)/u(1)
  w(3) = u(3)/u(1)
  w(4) = (gamma-one)*( u(4) - half*w(1)*(w(2)*w(2)+w(3)*w(3)) )

 end function u2w
!--------------------------------------------------------------------------------



!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!* Below, you find the following subroutines associated with LSQ gradients:
!*
!*  - compute_lsq_coeff_nc : Compute and store LSQ gradient coefficients at node
!*  - check_lsq_coeff_nc   : Check the computed LSQ coefficients
!*  - compute_gradient_nc  : Compute gradients at nodes
!*
!*  - lsq_gradients_nc     : Compute gradients (loop over neighbros) at a node
!*  - lsq_gradients2_nc    : Compute gradients (including neighbors of neighbors) at a node
!*  - lsq01_2x2_coeff_nc   : Compute linear LSQ coefficients at a node
!*  - lsq02_5x5_coeff2_nc  : Compute quadratic LSQ coefficients at all nodes
!*  - lsq_weight           : Compute LSQ weights
!*  - gewp_solve           : Gauss elimination to ivert a matrix
!
!********************************************************************************
!********************************************************************************

!********************************************************************************
!* This subroutine computes the coefficients for linear and quadratic LSQ
!* gradeints. 
!*
!* Note: Quadratic LSQ method is implemented in two steps as described in
!*       Nishikawa, JCP2014v273pp287-309 for details, which is available at 
!* http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2014v273pp287-309_preprint.pdf.
!*       This two-step method is useful in a paralell environment: it takes into
!*       account neighbors of neighbors without accessing neighbors of neighbors.
!*       So, it can be implemented only with edge-connected-neighbor information.
!*       In other words, each step is compact.
!*
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* This is Version 0 (July 2015).
!* This F90 code is written and made available for an educational purpose.
!* This file may be updated in future.
!*
!********************************************************************************
 subroutine compute_lsq_coeff_nc

 use edu2d_my_main_data , only : nnodes, node, nq
 use edu2d_my_allocation, only : my_alloc_p2_ptr, my_alloc_p2_matrix_ptr

 integer :: i, in, ell, ii, k

  write(*,*)
  write(*,*) "Constructing LSQ coefficients..."

! 1. Coefficients for the linear LSQ gradients

  write(*,*) "---(1) Constructing Linear LSQ coefficients..."

  nodes : do i = 1, nnodes

     call my_alloc_p2_ptr(node(i)%lsq2x2_cx,node(i)%nnghbrs)
     call my_alloc_p2_ptr(node(i)%lsq2x2_cy,node(i)%nnghbrs)
     call lsq01_2x2_coeff_nc(i)

  end do nodes

! 2. Coefficients for the quadratic LSQ gradients (two-step method)

  write(*,*) "---(2) Constructing Quadratic LSQ coefficients..."

  do i = 1, nnodes

        ii = 0 
     nghbr : do k = 1, node(i)%nnghbrs
        in = node(i)%nghbr(k)
      nghbr_nghbr : do ell = 1, node(in)%nnghbrs
        ii = ii + 1
      end do nghbr_nghbr
     end do nghbr

     call my_alloc_p2_ptr(node(i)%lsq5x5_cx, ii)
     call my_alloc_p2_ptr(node(i)%lsq5x5_cy, ii)
     call my_alloc_p2_ptr(node(i)%dx,node(i)%nnghbrs)
     call my_alloc_p2_ptr(node(i)%dy,node(i)%nnghbrs)
     call my_alloc_p2_matrix_ptr(node(i)%dw, nq,node(i)%nnghbrs)

  end do

     call lsq02_5x5_coeff2_nc

 end subroutine compute_lsq_coeff_nc

!********************************************************************************
!* This subroutine verifies the implementation of LSQ gradients.
!*
!* 1. Check if the linear LSQ gradients are exact for linear functions.
!* 2. Check if the quadratic LSQ gradients are exact for quadratic functions.
!*
!* Note: Here, we use only the first component of u=(u1,u2,u3), i.e., ivar=1.
!*
!********************************************************************************
 subroutine check_lsq_coeff_nc

 use edu2d_constants   , only : p2, one, two
 use edu2d_my_main_data, only : nnodes, node

 integer       :: i, ix, iy, ivar
 character(80) :: grad_type_temp
 real(p2)      :: error_max_wx, error_max_wy, x, y
 real(p2)      :: x_max_wx, y_max_wx, x_max_wy, y_max_wy, wx, wxe, wy, wye
 real(p2)      :: a0, a1, a2, a3, a4, a5

  ix = 1
  iy = 2

! We only use w(1) for this test.
  ivar = 1
 
!---------------------------------------------------------------------
! 1. Check linear LSQ gradients
!---------------------------------------------------------------------
  write(*,*)
  write(*,*) "---------------------------------------------------------"
  write(*,*) "---------------------------------------------------------"
  write(*,*) "- Checking Linear LSQ gradients..."

!  (1). Store a linear function in w(ivar) = x + 2*y.
!       So the exact gradient is grad(w(ivar)) = (1,2).

   write(*,*) "- Storing a linear function values..."
   do i = 1, nnodes
    x = node(i)%x
    y = node(i)%y
    node(i)%w(ivar) = one*x + two*y
   end do

!  (2). Compute the gradient by linear LSQ

   write(*,*) "- Computing linear LSQ gradients.."
   grad_type_temp = 'linear'
   call compute_gradient_nc(ivar,grad_type_temp)

!  (3). Compute the relative errors (L_infinity)

   write(*,*) "- Computing the relative errors (L_infinity).."
   error_max_wx = -one
   error_max_wy = -one
   do i = 1, nnodes
    error_max_wx = max( abs( node(i)%gradw(ivar,ix) - one )/one, error_max_wx )
    error_max_wy = max( abs( node(i)%gradw(ivar,iy) - two )/two, error_max_wy )
   end do

  write(*,*) " Max relative error in wx = ", error_max_wx
  write(*,*) " Max relative error in wy = ", error_max_wy
  write(*,*) "---------------------------------------------------------"
  write(*,*) "---------------------------------------------------------"


!---------------------------------------------------------------------
! 2. Check quadratic LSQ gradients
!---------------------------------------------------------------------
  write(*,*)
  write(*,*) "---------------------------------------------------------"
  write(*,*) "---------------------------------------------------------"
  write(*,*) "- Checking Quadratic LSQ gradients..."

!  (1). Store a quadratic function in w(ivar) = a0 + a1*x + a2*y + a3*x**2 + a4*x*y + a5*y**2
!       So the exact gradient is grad(w(ivar)) = (a1+2*a3*x+a4*y, a2+2*a5*y+a4*x)

   a0 =     21.21_p2
   a1 =     1.000_p2
   a2 = -   1.970_p2
   a3 =   280.400_p2
   a4 = -2129.710_p2
   a5 =   170.999_p2

   write(*,*) "- Storing a quadratic function values..."
   do i = 1, nnodes
    x = node(i)%x
    y = node(i)%y
    node(i)%w(ivar) = a0 + a1*x + a2*y + a3*x**2 + a4*x*y + a5*y**2
   end do

!  (2). Compute the gradient by linear LSQ

   write(*,*) "- Computing quadratic LSQ gradients.."
   grad_type_temp = 'quadratic2'
   call compute_gradient_nc(ivar,grad_type_temp)

!  (3). Compute the relative errors (L_infinity)

   write(*,*) "- Computing the relative errors (L_infinity).."
   error_max_wx = -one
   error_max_wy = -one
   do i = 1, nnodes
    x = node(i)%x
    y = node(i)%y

    if ( abs( node(i)%gradw(ivar,ix) - (a1+2.0_p2*a3*x+a4*y) )/(a1+2.0_p2*a3*x+a4*y) >  error_max_wx ) then
      wx  = node(i)%gradw(ivar,ix)
      wxe = a1+2.0_p2*a3*x+a4*y
      error_max_wx = abs( wx - wxe )/wxe
      x_max_wx = x
      y_max_wx = y
    endif

    if ( abs( node(i)%gradw(ivar,iy) - (a2+2.0_p2*a5*y+a4*x) )/(a2+2.0_p2*a5*y+a4*x) >  error_max_wy ) then
      wy  = node(i)%gradw(ivar,iy)
      wye = a2+2.0_p2*a5*y+a4*x
      error_max_wy = abs( wy - wye )/wye
      x_max_wy = x
      y_max_wy = y
    endif

   end do

  write(*,'(a,es20.3,a,2es12.5)') " Max relative error in wx = ", error_max_wx, " at (x,y) = ", x_max_wx, y_max_wx
  write(*,'(a,es20.10,a,es20.10)')  "   At this location, LSQ ux = ", wx, ": Exact ux = ", wxe
  write(*,'(a,es20.3,a,2es12.5)') " Max relative error in wy = ", error_max_wy, " at (x,y) = ", x_max_wy, y_max_wy
  write(*,'(a,es20.10,a,es20.10)')  "   At this location, LSQ uy = ", wy, ": Exact uy = ", wye
  write(*,*) "---------------------------------------------------------"
  write(*,*) "---------------------------------------------------------"
  write(*,*)


 end subroutine check_lsq_coeff_nc

!********************************************************************************
!* This subroutine computes gradients at nodes for the variable u(ivar),
!* where ivar = 1,2,3, ..., or nq.
!*
!* ------------------------------------------------------------------------------
!*  Input: node(:)%u(ivar)
!*
!* Output: node(i)%gradu(ivar,1:2) = ( du(ivar)/dx, du(ivar)/dy )
!* ------------------------------------------------------------------------------
!********************************************************************************
 subroutine compute_gradient_nc(ivar,grad_type)

 use edu2d_my_main_data, only : node, nnodes

 integer, intent(in) :: ivar

 integer       :: i, k, in
 character(80) :: grad_type

  if (trim(grad_type) == "none") return

!-------------------------------------------------
!  Two-step quadratic LSQ 5x5 system
!  Note: See Nishikawa, JCP2014v273pp287-309 for details, which is available at 
!        http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2014v273pp287-309_preprint.pdf.

   if (trim(grad_type) == "quadratic2") then

!  Perform Step 1 as below (before actually compute the gradient).

      do i = 1, nnodes

       nghbr0 : do k = 1, node(i)%nnghbrs
                   in      = node(i )%nghbr(k)
        node(i)%dx(k)      = node(in)%x       - node(i)%x
        node(i)%dy(k)      = node(in)%y       - node(i)%y
        node(i)%dw(ivar,k) = node(in)%w(ivar) - node(i)%w(ivar)
       end do nghbr0
      end do

   endif
!-------------------------------------------------

!------------------------------------------------------------
!------------------------------------------------------------
!-- Compute LSQ Gradients at all nodes.
!------------------------------------------------------------
!------------------------------------------------------------

  nodes : do i = 1, nnodes

  !-------------------------------------------------
  ! Linear LSQ 2x2 system
    if (trim(grad_type) == "linear") then

       call lsq_gradients_nc(i,ivar)

  !-------------------------------------------------
  ! Two-step quadratic LSQ 5x5 system
  !  Note: See Nishikawa, JCP2014v273pp287-309 for details, which is available at 
  !        http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2014v273pp287-309_preprint.pdf.
    elseif (trim(grad_type) == "quadratic2") then

      call lsq_gradients2_nc(i,ivar)

  !-------------------------------------------------
    else

     write(*,*) " Invalid input value -> ", trim(grad_type)
     stop

    endif
  !-------------------------------------------------

   end do nodes

 end subroutine compute_gradient_nc


!********************************************************************************
!* Compute the gradient, (wx,wy), for the variable u by Linear LSQ.
!*
!* ------------------------------------------------------------------------------
!*  Input:            inode = Node number at which the gradient is computed.
!*                     ivar =   Variable for which the gradient is computed.
!*          node(:)%w(ivar) = Solution at nearby nodes.
!*
!* Output:  node(inode)%gradu = gradient of the requested variable
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine lsq_gradients_nc(inode,ivar)

 use edu2d_my_main_data           , only : node
 use edu2d_constants              , only : p2, zero

 implicit none

 integer, intent(in) :: inode, ivar

!Local variables
 integer  :: in, inghbr
 integer  :: ix, iy
 real(p2) :: da, ax, ay

  ix = 1 
  iy = 2
  ax = zero
  ay = zero

!   Loop over neighbors

     do in = 1, node(inode)%nnghbrs
       inghbr = node(inode)%nghbr(in)

          da = node(inghbr)%w(ivar) - node(inode)%w(ivar)
 
      ax = ax + node(inode)%lsq2x2_cx(in)*da
      ay = ay + node(inode)%lsq2x2_cy(in)*da

     end do

      node(inode)%gradw(ivar,ix) = ax  !<-- dw(ivar)/dx
      node(inode)%gradw(ivar,iy) = ay  !<-- dw(ivar)/dy

 end subroutine lsq_gradients_nc
!--------------------------------------------------------------------------------


!********************************************************************************
!* Compute the gradient, (wx,wy), for the variable u by Quadratic LSQ.
!*
!* ------------------------------------------------------------------------------
!*  Input:            inode = Node number at which the gradient is computed.
!*                     ivar =   Variable for which the gradient is computed.
!*          node(:)%w(ivar) = Solution at nearby nodes.
!*
!* Output:  node(inode)%gradu = gradient of the requested variable
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine lsq_gradients2_nc(inode,ivar)

 use edu2d_my_main_data, only : node
 use edu2d_constants   , only : p2, zero

 implicit none

 integer, intent(in) :: inode, ivar

!Local variables
 integer  :: in
 integer  :: ix, iy, ii, ell, k
 real(p2) :: da, ax, ay

   ix = 1 
   iy = 2
   ax = zero
   ay = zero

!   Loop over neighbors

       ii = 0

     nghbr : do k = 1, node(inode)%nnghbrs
              in = node(inode)%nghbr(k)

      nghbr_nghbr : do ell = 1, node(in)%nnghbrs

       da = node(in)%w(ivar) - node(inode)%w(ivar) + node(in)%dw(ivar,ell)

       if ( node(in)%nghbr(ell) == inode ) then
        da = node(in)%w(ivar) - node(inode)%w(ivar)
       endif

       ii = ii + 1

       ax = ax + node(inode)%lsq5x5_cx(ii)*da
       ay = ay + node(inode)%lsq5x5_cy(ii)*da

      end do nghbr_nghbr

     end do nghbr

      node(inode)%gradw(ivar,ix) = ax  !<-- dw(ivar)/dx
      node(inode)%gradw(ivar,iy) = ay  !<-- dw(ivar)/dy

 end subroutine lsq_gradients2_nc
!--------------------------------------------------------------------------------



!********************************************************************************
!* --- LSQ Coefficients for 2x2 Linear Least-Squares Gradient Reconstruction ---
!*
!* ------------------------------------------------------------------------------
!*  Input:  inode = node number at which the gradient is computed.
!*
!* Output:  node(inode)%lsq2x2_cx(:)
!*          node(inode)%lsq2x2_cx(:)
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine lsq01_2x2_coeff_nc(inode)

 use edu2d_my_main_data           , only : node
 use edu2d_constants              , only : p2, zero

 implicit none

 integer, intent(in) :: inode
!Local variables
 real(p2) :: a(2,2), dx, dy, det, w2, w2dvar
 integer  :: k, inghbr, ix=1,iy=2
 real(p2), dimension(2,2) :: local_lsq_inverse

   a = zero

!  Loop over the neighbor nodes.
   do k = 1, node(inode)%nnghbrs
    inghbr = node(inode)%nghbr(k)

      dx = node(inghbr)%x - node(inode)%x
      dy = node(inghbr)%y - node(inode)%y

      w2 = lsq_weight(dx, dy)**2

      a(1,1) = a(1,1) + w2 * dx*dx
      a(1,2) = a(1,2) + w2 * dx*dy

      a(2,1) = a(2,1) + w2 * dx*dy
      a(2,2) = a(2,2) + w2 * dy*dy

   end do

    det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
    if (abs(det) < 1.0e-14) write(*,*) " Singular: LSQ det = ", det, " i=",inode

! OK, invert and store the inverse matrix:

     local_lsq_inverse(1,1) =  a(2,2)/det
     local_lsq_inverse(1,2) = -a(2,1)/det
     local_lsq_inverse(2,1) = -a(1,2)/det
     local_lsq_inverse(2,2) =  a(1,1)/det

!  Now compute the coefficients for neighbors.

     nghbr : do k = 1, node(inode)%nnghbrs
       inghbr = node(inode)%nghbr(k)

        dx = node(inghbr)%x - node(inode)%x
        dy = node(inghbr)%y - node(inode)%y

      w2dvar = lsq_weight(dx, dy)**2

      node(inode)%lsq2x2_cx(k)  = local_lsq_inverse(ix,1)*w2dvar*dx &
                                + local_lsq_inverse(ix,2)*w2dvar*dy

      node(inode)%lsq2x2_cy(k)  = local_lsq_inverse(iy,1)*w2dvar*dx &
                                + local_lsq_inverse(iy,2)*w2dvar*dy

     end do nghbr

 end subroutine lsq01_2x2_coeff_nc
!********************************************************************************
!*
!********************************************************************************


!********************************************************************************
!* --- LSQ Coefficients for 5x5 Quadratic Least-Squares Gradient Reconstruction ---
!*
!* Note: See Nishikawa, JCP2014v273pp287-309 for details, which is available at 
!*
!* http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2014v273pp287-309_preprint.pdf.
!*
!* ------------------------------------------------------------------------------
!*  Input:
!*
!* Output:  node(:)%lsq5x5_cx(ii)
!*          node(:)%lsq5x5_cx(ii)
!*
!* Note: This subroutine computes the LSQ coefficeints at all nodes.
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine lsq02_5x5_coeff2_nc

 use edu2d_my_main_data           , only : node, nnodes
 use edu2d_constants              , only : p2, zero, half

 implicit none

!Local variables
 real(p2) :: a(5,5), ainv(5,5), dx, dy
 real(p2) :: dummy1(5), dummy2(5)
 real(p2) :: w2
 integer  :: istat, ix=1, iy=2

 integer :: i, k, ell, in, ii

! Step 1

  node1 : do i = 1, nnodes
   nghbr : do k = 1, node(i)%nnghbrs

               in   = node(i)%nghbr(k)
    node(i)%dx(k)   = node(in)%x - node(i)%x
    node(i)%dy(k)   = node(in)%y - node(i)%y

   end do nghbr
  end do node1

! Step 2

  node2 : do i = 1, nnodes

     a = zero

!  Get dx, dy, and dw

   nghbr2 : do k = 1, node(i)%nnghbrs

              in = node(i)%nghbr(k)

    nghbr_nghbr : do ell = 1, node(in)%nnghbrs

     dx = node(in)%x - node(i)%x + node(in)%dx(ell)
     dy = node(in)%y - node(i)%y + node(in)%dy(ell)

     if ( node(in)%nghbr(ell) == i ) then

      dx = node(in)%x - node(i)%x
      dy = node(in)%y - node(i)%y

      if ( abs(dx) + abs(dy) < 1.0e-13_p2 ) then
       write(*,*) " Zero distance found at lsq02_5x5_coeff2_nc..."
       write(*,*) "    dx = ", dx
       write(*,*) "    dy = ", dy
       write(*,*) "- Centered node = ", i
       write(*,*) "          (x,y) = ", node(in)%x, node(in)%y
       write(*,*) "- Neighbor node = ", in
       write(*,*) "          (x,y) = ", node(in)%x, node(in)%y
       stop
      endif

     endif


        w2 = lsq_weight(dx, dy)**2

      a(1,1) = a(1,1) + w2 * dx         *dx
      a(1,2) = a(1,2) + w2 * dx         *dy
      a(1,3) = a(1,3) + w2 * dx         *dx*dx * half
      a(1,4) = a(1,4) + w2 * dx         *dx*dy
      a(1,5) = a(1,5) + w2 * dx         *dy*dy * half

!     a(2,1) = a(2,1) + w2 * dy         *dx
      a(2,2) = a(2,2) + w2 * dy         *dy
      a(2,3) = a(2,3) + w2 * dy         *dx*dx * half
      a(2,4) = a(2,4) + w2 * dy         *dx*dy
      a(2,5) = a(2,5) + w2 * dy         *dy*dy * half

!     a(3,1) = a(3,1) + w2 * half*dx*dx *dx
!     a(3,2) = a(3,2) + w2 * half*dx*dx *dy
      a(3,3) = a(3,3) + w2 * half*dx*dx *dx*dx * half
      a(3,4) = a(3,4) + w2 * half*dx*dx *dx*dy
      a(3,5) = a(3,5) + w2 * half*dx*dx *dy*dy * half

!     a(4,1) = a(4,1) + w2 *      dx*dy *dx
!     a(4,2) = a(4,2) + w2 *      dx*dy *dy
!     a(4,3) = a(4,3) + w2 *      dx*dy *dx*dx * half
      a(4,4) = a(4,4) + w2 *      dx*dy *dx*dy
      a(4,5) = a(4,5) + w2 *      dx*dy *dy*dy * half

!     a(5,1) = a(5,1) + w2 * half*dy*dy *dx
!     a(5,2) = a(5,2) + w2 * half*dy*dy *dy
!     a(5,3) = a(5,3) + w2 * half*dy*dy *dx*dx * half
!     a(5,4) = a(5,4) + w2 * half*dy*dy *dx*dy
      a(5,5) = a(5,5) + w2 * half*dy*dy *dy*dy * half

    end do nghbr_nghbr

   end do nghbr2

!   Fill symmetric part

      a(2,1) = a(1,2);
      a(3,1) = a(1,3);  a(3,2) = a(2,3);
      a(4,1) = a(1,4);  a(4,2) = a(2,4); a(4,3) = a(3,4);
      a(5,1) = a(1,5);  a(5,2) = a(2,5); a(5,3) = a(3,5); a(5,4) = a(4,5);

!   Invert the matrix

     dummy1 = zero
     dummy2 = zero
     call gewp_solve(a,dummy1,dummy2,ainv,istat, 5);

     if (istat/=0) write(*,*) "Problem in solving the linear system!: Quadratic_LSJ_Matrix"

!  Now compute the coefficients for neighbors.

      ii = 0

     nghbr3 : do k = 1, node(i)%nnghbrs
                in = node(i)%nghbr(k)

      nghbr_nghbr2 : do ell = 1, node(in)%nnghbrs

       dx = node(in)%x - node(i)%x + node(in)%dx(ell)
       dy = node(in)%y - node(i)%y + node(in)%dy(ell)

       if ( node(in)%nghbr(ell) == i ) then
        dx = node(in)%x - node(i)%x
        dy = node(in)%y - node(i)%y
       endif

       ii = ii + 1
       
       w2 = lsq_weight(dx, dy)**2

 !  Multiply the inverse LSQ matrix to get the coefficients: cx(:) and cy(:):

      node(i)%lsq5x5_cx(ii)  =  ainv(ix,1)*w2*dx             &
                              + ainv(ix,2)*w2*dy             &
                              + ainv(ix,3)*w2*dx*dx * half   &
                              + ainv(ix,4)*w2*dx*dy          &
                              + ainv(ix,5)*w2*dy*dy * half

      node(i)%lsq5x5_cy(ii)  =  ainv(iy,1)*w2*dx             &
                              + ainv(iy,2)*w2*dy             &
                              + ainv(iy,3)*w2*dx*dx * half   &
                              + ainv(iy,4)*w2*dx*dy          &
                              + ainv(iy,5)*w2*dy*dy * half
      end do nghbr_nghbr2

     end do nghbr3

  end do node2

 end subroutine lsq02_5x5_coeff2_nc
!********************************************************************************
!*
!********************************************************************************


!****************************************************************************
!* Compute the LSQ weight
!*
!* Note: The weight computed here is the square of the actual LSQ weight.
!*****************************************************************************
 function lsq_weight(dx, dy)

 use edu2d_constants   , only : p2, one
 use edu2d_my_main_data, only : gradient_weight, gradient_weight_p

 implicit none

!Input
 real(p2), intent(in) :: dx, dy
!Output
 real(p2)             :: lsq_weight
!Local
 real(p2)             :: distance

  if     (trim(gradient_weight) == "none"            ) then

        lsq_weight = one

  elseif (trim(gradient_weight) == "inverse_distance") then

          distance = sqrt(dx*dx + dy*dy)
        lsq_weight = one / distance**gradient_weight_p

  endif

 end function lsq_weight

!****************************************************************************
!* ------------------ GAUSS ELIMINATION WITH PIVOTING ---------------------
!*
!*  This computes the inverse of an (Nm)x(Nm) matrix "ai".
!*
!*  IN :       ai = an (Nm)x(Nm) matrix whoise inverse is sought.
!*             bi = right hand side of a linear system: ai*x=bi.
!*             nm = the size of the matrix "ai"
!*
!* OUT :  inverse = the inverse of "ai".
!*            sol = solution to the linear system, ai*x=bi
!*       idetstat = 0 -> inverse successfully computed
!*                  1 -> THE INVERSE DOES NOT EXIST (det=0).
!*                  2 -> No unique solutions exist.
!*
!* Katate Masatsuka, April 2015. http://www.cfdbooks.com
!*****************************************************************************
 subroutine gewp_solve(ai,bi,sol,inverse,idetstat,nm)

  use edu2d_constants   , only : p2, zero, one

  implicit none

! Input
  real(p2), intent( in) :: ai(:,:),bi(:)
  integer , intent( in) :: nm

! Output
  real(p2), intent(out) :: sol(:),inverse(nm,nm)
  integer , intent(out) :: idetstat

! Local variables
  real(p2) :: a(nm,nm+1),x(nm)
  integer  :: i,j,k,pp,nrow(nm),m

 do m=1,nm
!*****************************************************************************
!*****************************************************************************
       do j=1,nm
        do i=1,nm
          a(i,j) = ai(i,j)
        end do
       end do
       do k=1,nm
          a(k,nm+1)=zero; nrow(k)=k
       end do
          a(m,nm+1)=one
!*****************************************************************************
       do j=1,nm-1
!*****************************************************************************
!* find smallest pp for a(pp,j) is maximum in jth column.
!***************************************************************************** 
      call findmax(nm,j,pp,a,nrow)
!*****************************************************************************
!* if a(nrow(p),j) is zero, there's no unique solutions      
!*****************************************************************************
      if ( abs(a(nrow(pp),j) - zero) < 1.0e-15 ) then
       write(6,*) 'the inverse does not exist.'
        idetstat = 1
        return
      else
      endif
!*****************************************************************************
!* if the max is not a diagonal element, switch those rows       
!*****************************************************************************
      if (nrow(pp) .ne. nrow(j)) then
      call switch(nm,j,pp,nrow)
      else
      endif  
!*****************************************************************************
!* eliminate all the entries below the diagonal one
!***************************************************************************** 
      call eliminate_below(nm,j,a,nrow)

      end do
!*****************************************************************************
!* check if a(nrow(n),n)=0.0 .
!*****************************************************************************
      if ( abs(a(nrow(nm),nm) - zero) < 1.0e-15 ) then
        write(6,*) 'no unique solution exists!';  idetstat = 2
        return 
      else
      endif
!*****************************************************************************
!* backsubstitution!
!*****************************************************************************
      call backsub(nm,x,a,nrow)
!*****************************************************************************
!* store the solutions, you know they are inverse(i,m) i=1...
!*****************************************************************************
      do i=1,nm
         inverse(i,m)=x(i)
      end do
!*****************************************************************************
 end do
!*****************************************************************************
!* solve
!*****************************************************************************
    do i=1,nm; sol(i)=zero;
     do j=1,nm
       sol(i) = sol(i) + inverse(i,j)*bi(j)
     end do
    end do

    idetstat = 0;
    return

 end subroutine gewp_solve


!Four subroutines below are used in gewp above.


!*****************************************************************************
!* Find maximum element in jth column 
!***************************************************************************** 
 subroutine findmax(nm,j,pp,a,nrow)

  use edu2d_constants   , only : p2

  implicit none

! Input
  integer , intent( in) :: nm
  real(p2), intent( in) :: a(nm,nm+1)
  integer , intent( in) :: j,nrow(nm)

! Output
  integer , intent(out) :: pp

! Local variables
  real(p2) :: max
  integer :: i

            max=abs(a(nrow(j),j)); pp=j
           do i=j+1,nm
             if (max < abs(a(nrow(i),j))) then
                  pp=i; max=abs(a(nrow(i),j))
             endif
           end do

  return

 end subroutine findmax

!*****************************************************************************
!* Switch rows       
!*****************************************************************************
 subroutine switch(nm,j,pp,nrow)

 implicit none

! Input
  integer, intent(   in) :: nm,j,pp

! Output
  integer, intent(inout) :: nrow(nm)

! Local
  integer :: ncopy

      if (nrow(pp).ne.nrow(j)) then
         ncopy=nrow(j)
         nrow(j)=nrow(pp)
         nrow(pp)=ncopy
      endif  

  return

 end subroutine switch

!*****************************************************************************
!* Eliminate all the entries below the diagonal one
!* (give me j, the column you are working on now)
!***************************************************************************** 
 subroutine eliminate_below(nm,j,a,nrow)

  use edu2d_constants   , only : p2, zero

  implicit none
  
! Input
  integer , intent(   in) :: nm
  integer , intent(   in) :: j,nrow(nm)

! Output
  real(p2), intent(inout) :: a(nm,nm+1)

! Local
  real(p2) :: m
  integer  :: k,i

      do i=j+1,nm
        m=a(nrow(i),j)/a(nrow(j),j)
        a(nrow(i),j)=zero
          do k=j+1,nm+1
            a(nrow(i),k)=a(nrow(i),k)-m*a(nrow(j),k)
          end do
      end do

  return

 end subroutine eliminate_below

!*****************************************************************************
!* Backsubstitution!
!*****************************************************************************
 subroutine backsub(nm,x,a,nrow)

  use edu2d_constants   , only : p2, zero

  implicit none

! Input
  integer , intent( in) :: nm
  real(p2), intent( in) :: a(nm,nm+1)
  integer , intent( in) :: nrow(nm)

! Output
  real(p2), intent(out) :: x(nm)

! Local
  real(p2) :: sum
  integer :: i,k

      x(nm)=a(nrow(nm),nm+1)/a(nrow(nm),nm)
      do i=nm-1,1,-1
         sum=zero
           do k=i+1,nm
              sum=sum+a(nrow(i),k)*x(k)
           end do
      x(i)=(a(nrow(i),nm+1)-sum)/a(nrow(i),i)
      end do

  return

 end subroutine backsub
!*********************************************************************

 end module edu2d_euler_implct_solver
