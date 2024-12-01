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
 use edu2d_euler_linear_solve, only : gs_sequential, gs_sequential2 ! GS relaxaton for linear system
 use residual, only : compute_residual_ncfv

 implicit none

!Local variables
 real(p2), dimension(nq,3)             :: res_norm      !Residual norms(L1,L2,Linf)
 real(p2), dimension(:,:), allocatable :: u0            !Saved solution
 integer                               :: i_iteration   !Iteration counter
 real(p2)                              :: s, exp_factor !CFL ramping variables
 integer                               :: i
 real(p2)                              :: time_begin    !Starting time
 real(p2)                              :: time_end      !End time
 integer                               :: sweeps_actual
 real(p2)                              :: roc

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
    !call gs_sequential           ! -> This relaxes the linear system and compute du.
    call gs_sequential2(sweeps_actual,roc) ! -> This relaxes the linear system and compute du.

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
 end module edu2d_euler_implct_solver
