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

 use edu2d_my_main_data, only : CFL, nnodes, node, tolerance, CFL1, CFL2, CFL_ramp_steps, &
                                iteration_method, max_iterations, CFLexp, nq, i_iteration, &
								max_projection_gcr

!Implicit method uses subroutines below to construct Jacobian matrix
!and relax the linear system.
 use edu2d_euler_jacobian    , only : construct_jacobian_ncfv ! Construct Jacobian
 use edu2d_euler_linear_solve, only : gs_sequential, gs_sequential2 ! GS relaxaton for linear system
 use gcr_solver, only : jfnk_gcr_solver
 use residual,   only : compute_residual_ncfv

 implicit none

!Local variables
 real(p2), dimension(nq) :: res_norm1_initial
 real(p2), dimension(nq) :: res_norm1_pre

 real(p2), dimension(nq,3)             :: res_norm      !Residual norms(L1,L2,Linf)
 real(p2), dimension(:,:), allocatable :: u0            !Saved solution
 real(p2)                              :: s, exp_factor !CFL ramping variables
 integer                               :: i
 real(p2)                              :: time_begin    !Starting time
 real(p2)                              :: time_end      !End time
 integer                               :: sweeps_actual, gcr_actual
 integer                               :: total_steps
 real(p2)                              :: roc, CFL_previous

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

  iteration : do i_iteration = 0, max_iterations-1

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
    if (i_iteration==0) then
      write(*,'(a94)') "Density    X-momentum  Y-momentum   Energy"
    endif

    if (mod(i_iteration,1)==0) write(*,'(a4,es12.2,a13,i6,a12,4es12.2)') &
       "CFL=", CFLexp, "Iteration=", i_iteration, " L1(res)=",res_norm(:,1)

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
       s = real( i_iteration, p2) / real( CFL_ramp_steps, p2)
     CFL = CFL1 + (CFL2-CFL1)*( one - exp(-s*exp_factor) )/ ( one - exp(-exp_factor) )

    endif

!   Compute Residual: Right hand side of [V/dt+dR/dU]*dU = -Residual
    call compute_residual_ncfv   ! -> This computes node(:)%res.

!   Compute residual norms (undivided residual)
    call residual_norm(res_norm) ! -> This computes res_norm.

!   Display the residual norm.
    if (i_iteration==0) then
      write(*,'(a101)') "Density    X-momentum  Y-momentum   Energy"
    endif

    if (mod(i_iteration,1)==0) write(*,'(a4,f20.2,a13,i6,a12,4es12.2)') &
       "CFL=", CFL, "iteration=", i_iteration, " L1(res)=",res_norm(:,1)

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
!*********************************************************************
!*********************************************************************
! Method 3: Implicit JFNK-GCR solver with DC preconditioner
!	
  elseif (trim(iteration_method) == "implicit_gcr") then
  
    CFL_previous = CFL

!   Exponential CFL ramping between CFL1 and CFL2 over 'CFL_ramp_steps' iterations.
!   Note: This is by no means an optimal ramping strategy.
    if (i_iteration <= CFL_ramp_steps) then

     exp_factor = 0.5_p2
       s = real( i_iteration, p2) / real( CFL_ramp_steps, p2)
     CFL = CFL1 + (CFL2-CFL1)*( one - exp(-s*exp_factor) )/ ( one - exp(-exp_factor) )

    endif

!   Compute Residual: Right hand side of [V/dt+dR/dU]*dU = -Residual
    call compute_residual_ncfv   ! -> This computes node(:)%res.

!   Compute residual norms (undivided residual)
    call residual_norm(res_norm) ! -> This computes res_norm.
	
!   Set the initial residual after the first iteration.
    if (i_iteration==1) then
      res_norm1_initial = res_norm(:,1)
      res_norm1_pre     = res_norm(:,1)
    endif

!   Display the residual norm.
    if (i_iteration==0) write(*,'(a101)') "Density    X-momentum  Y-momentum   Energy"
	if (i_iteration>=0) then
     call report_res_norm(CFL_previous,res_norm(:,:),res_norm1_pre, sweeps_actual,gcr_actual,roc)

     if (i_iteration==1) write(*,*)
     if (i_iteration==1) write(*,'(a)') "    Note: We take these residuals (after the 1st iteration)"
     if (i_iteration==1) write(*,'(a)') "          as the initial residual to measure the reduction."
     if (i_iteration==1) write(*,*)

     write(*,'(a86)') "   ----------------------------------------------------------------------------------"
    endif
	
!   Stop if the L1 residual norm drops below the specified tolerance or if we reach the max iteration.
!   Note: Machine zero residual level is very hard to define. Here, I just set it to be 1.0e-14.
!         It may be larger or smaller for different problems.
    if (i_iteration >= 1) then

     if (     maxval(res_norm(:,1)/res_norm1_initial(:)) < tolerance & !Tolerance met
         .or. maxval(res_norm(:,1)) < 1.0e-14                        & !Machine zero(?)
                                                                                 ) then
       total_steps = i_iteration
       write(*,*)
       write(*,'(a15,i10,a12,4es12.2,a10)') "    steps=", total_steps, &
                                       " L1(res)=",res_norm(:,1), " Converged"
       write(*,*)
       exit iteration
     endif

    endif
	
	!   Proceed to the next iteration
    res_norm1_pre = res_norm(:,1)

!   Compute the local time step (the local time step is used by DC in pseudo time term)
    call compute_time_step       ! -> This computes node(:)%dt.

!   Construct the residual Jacobian matrix (based on 1st-order scheme: Defect Correction)
!   which is the matrix on the left hand side of [V/dt+dR/dU]*dU = -Residual.
!   Note: This sbroutine can be found in edu2d_euler_jacobian_v0.f90
    call construct_jacobian_ncfv ! -> This computes jac(:)%diag and jac(:)%off.

!   Compute du (Correction) by Defect-correction iteration or Newton-Krylov.

    !------------------------------------------------
    !(1)Defect correction: Relax the linear system by Gauss-Seidel to get du (Correction)
     if (max_projection_gcr == 0) then

      call gs_sequential2(sweeps_actual,roc)

    !----------------------------------------------------
    !(2)GCR with defect-correction as a preconditioner
     elseif (max_projection_gcr > 0) then

      !For GCR, sweeps_actual=actual_projections, cr_gs=l1norm_ratio
       call jfnk_gcr_solver(gcr_actual,sweeps_actual,roc)

     endif
    !------------------------------------------------

!   Update the solution: u_new = u_old + du
    call update_solution_du      ! -> This updates node(:)%u.
	
	total_steps = i_iteration

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

  if (i_iteration == max_iterations-1) then
   write(*,*) " Not converged... Sorry..."
   write(*,*) "   Max iterations reached... max_iterations=", max_iterations
   write(*,*) "   Increase max_iterations, and try again."
  endif

  write(*,*) " Converged."
  write(*,*) " Final iteration      =", i_iteration
  write(*,*)
  write(*,*) "Finished the Euler solver... Bye!"

 end subroutine euler_solver_main
!--------------------------------------------------------------------------------

!********************************************************************************
!* This subroutine print out in the screen the residual history and associated
!* data.
!*
!********************************************************************************
 subroutine report_res_norm(CFL,res_norm,res_norm1_pre,sweeps,gcr,roc)

 use edu2d_constants   , only : p2
 use edu2d_my_main_data, only : i_iteration, nq, max_projection_gcr

 implicit none
 real(p2)                 , intent(in) :: CFL, roc
 real(p2), dimension(nq,3), intent(in) :: res_norm
 real(p2), dimension(nq)  , intent(in) :: res_norm1_pre

 integer                               :: k, sweeps, gcr
 real(p2), dimension(nq)               :: cr

 !Below, we print out only the first 3 components of res_norm and cr
 !because the equations we consider in this code has 3 equations.

   write(*,'(a10,es9.2, a9,i10,a12,4es12.2)') &
           "CFL=", CFL, "steps=", i_iteration, " L1(res)=",res_norm(:,1)

 !Compute the convergence rate:
 ! The ratio of the current residual norm to the previous one for each equation.
 ! So, it should be less than 1.0 if converging.

   do k = 1, nq

    if (res_norm1_pre(k) > 1.0e-15) then
     cr(k) = res_norm(k,1)/res_norm1_pre(k)
    else
     cr(k) = res_norm(k,1) !To avoid zero division.
    endif

   end do

 !Print convergence rate ( = the ratio of the current residual norm to the previous one for each
 !                             equation. So, it should be less than 1.0 if converging.)
  if (i_iteration > 1) write(*,'(33x,a15,4f12.4)') "    c.r.: ", cr(1:nq)

 !Print the number of sweeps/projections and the convergence rate.

   !(1)Defect-Correction case (The case when max_projection_gcr  =0)
   if (max_projection_gcr == 0) then

    if (i_iteration > 0) write(*,'(27x,a14,i5,a,es8.2)') "GS(sweeps:cr)=", sweeps,":",roc

   !(2)Newton-Krylov(GCR) case (The case when max_projection_gcr > 0)
   elseif (max_projection_gcr > 0) then

    if (i_iteration > 0) write(*,'(27x,a28,i5,a,i12,a,es8.2)') &
                         "GCR(projections:sweeps:cr)=", gcr,":",sweeps,":",roc
   else

    ! max_projection_gcr < 0 has no meaning... Stop.

    write(*,*) " Invalid value: max_projection_gcr = ", max_projection_gcr
    write(*,*) " Set zero or positive integer, and try again. Stop..."
    stop

   endif

 end subroutine report_res_norm
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
 use vector_operations,  only : u2w, w2u

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
 use vector_operations,  only : u2w, w2u

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
 end module edu2d_euler_implct_solver
