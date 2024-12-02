!--------------------------------------------------------------------------------


!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!* This subroutine performs one iteration of a Jabobian-Free Newton-Krylov method
!* based on GCR with a variable preconditioner.
!*
!* Note: This is a nonlinear solver, equivalent to Newton's mehtod.
!*       But the Jacobian is not computed nor stored.
!*       This is a Jacobian-Free Newton-Krylov method.
!*       The linear solver is a Krylov method (GCR).
!*       It requries a good preconditioner, for which we employ
!*       the implicit DC solver. This is a variable preconditioner
!*       because the preconditioning is not exactly the same every time.
!*       (A popular preconditioner is ILU, which is exactly the same every time.)
!*
!*
!* We solve the following linear system:
!*
!*    dRes/dU*dU = -Res    (A*x=b)
!*
!* The left hand side can be formed without computing and storing dRes/dU:
!*
!*  dRes/dU*dU ~ [ Res(U+epsilon*dU/|dU|) - Res(U) ] / [ epsilon*|dU| ]
!*
!* What is required in the GCR is the computation of dRes/dU*p, and so
!*
!*  dRes/dU*p ~ [ Res(U+epsilon*p/|p|) - Res(U) ] / [ epsilon*|p| ]
!*
!* ------------------------------------------------------------------------------
!*  Input:  node(:)%u   - Solution
!*          node(:)%res - Residual
!*
!*
!* Output:  node(:)%du  - Correction to be applied to the solution:
!*                        node(:)%u = node(:)%u + node(:)%du
!* ------------------------------------------------------------------------------
!*
!* Written by Katate Masatsuka (http://www.cfdbooks.com)
!********************************************************************************
module gcr_solver

private

public :: jfnk_gcr_solver

contains

 subroutine jfnk_gcr_solver(actual_projections,actual_sweeps,l2norm_ratio)

 use edu2d_constants   , only : p2, zero
 use edu2d_my_main_data, only : nq, nnodes, my_eps, sweeps_actual_gcr
 use input_parameter,    only : max_projection_gcr, tolerance_gcr
 use frechet_derivative, only : compute_Ap
 use precond, only : preconditioner
 use vector_operations, only : vector_res, elltwo_norm, length, vector2du

 integer , intent(out) :: actual_projections, actual_sweeps
 real(p2), intent(out) :: l2norm_ratio
 
 !Max dimension is set here (65 is too many, and so large enough!).
 integer , parameter                            :: max_projections=65
 real(p2), dimension(nnodes*nq)                 :: x, r, r0, a, b, Apip1
 real(p2), dimension(nnodes*nq,max_projections) :: p, Ap 

 real(p2) :: alpha, beta, l2norm_initial, length_p, length_r
 integer  :: i, k, ndim, istat_preconditioner

              alpha = 1.0_p2
  sweeps_actual_gcr = 0

  write(2000,*) ">>>"
  write(2000,*) ">>> GCR begins"

  if (max_projection_gcr > max_projections) then 
   write(*,*) "Too large max_projection_gcr..., it should be <", max_projections
   write(*,*) "Try again with a smaller number, or modify the code. Stop."
   stop
  endif

  ndim = nnodes*nq

! Initialization

   x = zero
   r = zero
   p = zero

!  Initial residual is the nonlinear residual since x=0 (r = b-A*x = b = -Res).

   r0 = vector_res(nnodes,nq) ! This is b=-Res.
   l2norm_initial = elltwo_norm(r0,ndim)
   write(2000,*) "i = ", 0, " L1(residual) = ", l2norm_initial

    ! p_{i+1} = A_approx^{-1}*r, where A_approx is an approximate Jacobian.
    call preconditioner(r0,ndim, p(:,1),istat_preconditioner)
    Ap(:,1) = compute_Ap(p(:,1),r0,ndim) !Ap by the Frechet derivative

     r  = r0

!---------------------------------------------------------------------
! Beginning of GCR Projection
!---------------------------------------------------------------------
   i = 0
  gcr_projection : do
   i = i + 1

    actual_projections = i

!   Here, for i > 1, alpha = beta = -(Ap_{i+1},Ap_{k})/(Ap_{k},Ap_{k})
    if (abs(alpha) < my_eps) then
     write(2000,*) ">>>>> Exit GCR: GCR Stall - Residual stalls with beta=", alpha
      x(:) = p(:,1) ! Cancel the update and exit: Revert to the defect correction.
     exit gcr_projection
    endif

!-----------------------------------------
!  Update x and r:
!
!      alpha = (r,Ap_i)/(Ap_i,Ap_i)  <- inner products.
!          x = x + alpha*p_i
!          r = r - alpha*Ap_i
!
!  Computation of alpha is implemented with normalized vectors as below:

    length_r = length(r,ndim)
    length_p = length(Ap(:,i),ndim)
           a =       r/max(my_eps,length_r)
           b = Ap(:,i)/max(my_eps,length_p)
       alpha = dot_product(a,b)*(length_r/length_p)
        x(:) = x(:) + alpha* p(:,i)
        r(:) = r(:) - alpha*Ap(:,i)

!-----------------------------------------
!  Check the convergence. If not converged, construct a new search direction.

     l2norm_ratio = elltwo_norm(r,ndim) / l2norm_initial
      write(2000,*) " i = ", i, " L1(residual) = ", l2norm_initial, " roc = ", l2norm_ratio

     if (l2norm_ratio < tolerance_gcr) then
      write(2000,*) ">>>>> Exit GCR: Converged..."
      exit gcr_projection
     endif

     if (l2norm_ratio > 100.0_p2) then
      write(2000,*) ">>>>> Exit GCR: Diverged over two orders of magnitude..."
      exit gcr_projection
     endif

     if (i == max_projection_gcr) then
      write(2000,*) ">>>>> Exit GCR: max_projection_gcr reached... max_projection_gcr = ", max_projection_gcr
      exit gcr_projection
     endif

     if (i == max_projections) then
      write(2000,*) ">>>>> Exit GCR: max_projections reached... max_projections = ", max_projections
      write(2000,*) ">>>>> This is the maximum projections specified inside the code..."
      exit gcr_projection
     endif

     if (l2norm_ratio > 0.999_p2) then
      write(2000,*) ">>>>> Exit GCR: l2norm_ratio = 1.0 (stall)... Revert to the defect correction."
      x(:) = p(:,1) ! Cancel the update and exit: Revert to the defect correction.
      exit gcr_projection
     endif

     if (l2norm_ratio > 0.99_p2 .and. i == 5) then
      write(2000,*) ">>>>> Exit GCR: l2norm_ratio = 0.99 after 5 projections. Revert to the defect correction."
      x(:) = p(:,1) ! Cancel the update and exit: Revert to the defect correction.
      exit gcr_projection
     endif

     if (istat_preconditioner > 0) then
      write(2000,*) ">>>>> Exit GCR: Preconditioner diverged."
      exit gcr_projection
     endif

!-----------------------------------------
!  Construct the next search direction.
!
!    p_{i+1} = A_approx^{-1}*r
!   Ap_{i+1} = [ Res(U+epsilon*p_{i+1}/|p_{i+1}|) - Res(U) ] / [ epsilon*|p_{i+1}| ]

    call preconditioner(r,ndim, p(:,i+1),istat_preconditioner)
    Ap(:,i+1) = compute_Ap(p(:,i+1),r0,ndim) !Ap by the Frechet derivative
        Apip1 = Ap(:,i+1)
    do k = 1, i

!      beta = (Ap_{i+1},Ap_{k})/(Ap_{k},Ap_{k})
!   p_{i+1} = p_{i+1} - beta*p_{k}
!
!   Computation of beta is done with normalized Ap vectors:

       length_p = length(Ap(:,k),ndim)
              a = Apip1  /max(my_eps,length_p)
              b = Ap(:,k)/max(my_eps,length_p)
           beta = dot_product(a,b)
      Ap(:,i+1) = Ap(:,i+1) - Ap(:,k)*beta !dot_product(a,b)
       p(:,i+1) =  p(:,i+1) -  p(:,k)*beta !dot_product(a,b)

    end do

  end do gcr_projection
!---------------------------------------------------------------------
! End of GCR Projection
!---------------------------------------------------------------------

  write(2000,*) " GCR ",i," ends.", " # of p vectors = ", actual_projections

! Store the solution x(:) in node(:)%du
  call vector2du(x,nnodes,nq)

  actual_sweeps = sweeps_actual_gcr

  write(2000,*) ">>> GCR ends."
  write(2000,*) ">>>"

 end subroutine jfnk_gcr_solver
!********************************************************************************
end module gcr_solver