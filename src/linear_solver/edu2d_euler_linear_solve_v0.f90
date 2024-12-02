!********************************************************************************
!* Educationally-Designed Unstructured 2D (EDU2D) Code
!*
!*
!*      This module belongs to the inviscid version: EDU2D-Euler-Implicit
!*
!*
!*
!* This module contains all subroutines needed to relax the linear system.
!*
!* Contained subroutines/functions:
!* -----------------------------------------------------------------------------
!* gs_sequential: Sequential Gauss Seidel relaxation
!*    gewp_solve: Gauss elimination with pivoting (to invert the diagonal block)
!* -----------------------------------------------------------------------------
!*
!* Written by Katate Masatsuka (http://www.cfdbooks.com)
!********************************************************************************
 module edu2d_euler_linear_solve

 private

 public :: gs_sequential
 public :: gs_sequential2
 public :: sgs_sequential

 contains

!********************************************************************************
!* This subroutine relaxes the linear system by Sequential Gauss-Seidel
!*
!* ------------------------------------------------------------------------------
!*  Input:         jac = Jacobian matirx (left hand side)
!*         node(:)%res = Residual (right hand side)
!*              sweeps = The number of GS sweeps to be performed
!*
!* Output: node(:)%du  = Correction
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine gs_sequential

 use edu2d_constants   , only : p2, zero
 use edu2d_my_main_data, only : nnodes, node, jac, sweeps, nq
 use gaussian_elimination, only : gewp_solve

 implicit none

!Local variables
 real(p2), dimension(nq,nq) :: inverse       ! Inverse of diagonal block
 real(p2), dimension(nq)    :: b, x
 integer                    :: i, k, idestat

 integer  :: j, ii, kth_nghbr

!--------------------------------------------------------------------------------
! 1. Initialize the correction

  nodes1 : do i = 1, nnodes
   node(i)%du = zero
  end do nodes1

!--------------------------------------------------------------------------------
! 2. Invert the diagonal block and store it (replace jac(i)%diag by its inverse)

   nodes2 : do i = 1, nnodes

!   b is just a dummy variable here.

     b = zero

!   Invert the diagonal block at node i by Gauss elimination with pivoting.
!
!   Note: gewp_solve() actually solves a linear system, Ax=b, but here
!         we use it only to obtain the inverse of A. So, b is a dummy.
!   Note: Gauss Elimination is not the only way. You can try LU decomposition,
!         iterative methods, etc.
!   Note: "x" is a dummy, not used.

!             Ax=b ->          A  b     x   A^{-1}
    call gewp_solve( jac(i)%diag, b, x, inverse, idestat, nq )

!   Report errors if it fails.

    if (idestat/=0) then
     write(*,*) " Error in inverting the diagonal block... Stop"
     do k = 1, nq
      write(*,'(4(es25.15))') ( jac(i)%diag(k,j), j=1,nq)
     end do
     stop
    endif

!   Save the inverse (replace "diag" by its inverse)

    jac(i)%diag = inverse

   end do nodes2

!--------------------------------------------------------------------------------
! 3. Linear Relaxation (Sweep)

  relax : do ii = 1, sweeps

!----------------------------------------------------
!    Sequential Gauss-Seidel Relaxation(sweep)

   nodes3 : do i = 1, nnodes

!    Form the right hand side of GS relaxation: [ sum( off_diagonal_block*du ) - residual ]

     b = node(i)%res ! Remember that residual has already been given the minus sign.

     neighbors_of_i : do k = 1, node(i)%nnghbrs

      kth_nghbr = node(i)%nghbr(k)
      b = b - matmul(jac(i)%off(:,:,k), node(kth_nghbr)%du) ! Matrix-vector multiplication

     end do neighbors_of_i

!    Update du by the GS relaxation
!    Note: Remember that diag is now the inverse of the diagonal block.

     node(i)%du = matmul( jac(i)%diag, b ) ! Matrix-vector multiplication

   end do nodes3

!    End of Sequential Gauss-Seidel Relaxation(sweep)
!----------------------------------------------------

  end do relax

!    End of Linear Relaxation (Sweep)
!--------------------------------------------------------------------------------

 end subroutine gs_sequential
!--------------------------------------------------------------------------------



!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!* This subroutine relaxes the linear system, Jac*du = -Res,
!* by the sequential Gauss-Seidel relaxation scheme.
!*
!* ------------------------------------------------------------------------------
!*  Input:         jac = Jacobian matirx (left hand side)
!*         node(:)%res = Residual (right hand side)
!*              sweeps = The number of GS sweeps to be performed
!*    tolerance_linear = Tolerance on the linear residual
!*
!*
!* Output: node(:)%du  = Correction
!*       sweeps_actual = Number of sweeps performed
!*                 roc = Rate of convergence (the ratio of the final linear residual
!*                       to the initial linear residual), showing how much it is
!*                       reduced. Tolerance is typically 1 order of magnitude.
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine gs_sequential2(sweeps_actual,roc)

 use edu2d_constants   , only : p2, zero
 use edu2d_my_main_data, only : nnodes, node, jac, sweeps, nq, &
                                tolerance_linear, i_iteration

 implicit none

 integer , intent(out)      :: sweeps_actual
 real(p2), intent(out)      :: roc

!Local variables
 real(p2), dimension(nq,nq) :: inverse
 real(p2), dimension(nq)    :: b, x
 integer                    :: i, k, idestat

 integer                    :: ii, kth_nghbr
 real(p2), dimension(nq,3)  :: linear_res_norm !Residual norms(L1,L2,Linf) for the linear system
 real(p2), dimension(nq,3)  :: linear_res_norm_initial
 real(p2), dimension(nq,3)  :: linear_res_norm_previous

 real(p2)                   :: omega           !Relaxation parameter

            omega = 1.0_p2

 write(1000,*) "DC Iteration = ", i_iteration

 sweeps_actual = 0

!--------------------------------------------------------------------------------
! 1. Initialize the correction

  nodes1 : do i = 1, nnodes
   node(i)%du = zero
  end do nodes1
!--------------------------------------------------------------------------------
! 2. Linear Relaxation (Sweep)

  relax : do ii = 1, sweeps

!----------------------------------------------------
!    Sequential Gauss-Seidel Relaxation (sweep)
   linear_res_norm(:,1) = zero

   nodes3 : do i = 1, nnodes

!   Form the right hand side of GS: [ sum( off_diagonal_block*du ) - residual ]

     b = node(i)%res ! Residual has already been given the minus sign.

    do k = 1, node(i)%nnghbrs

     kth_nghbr = node(i)%nghbr(k)
     b = b - matmul(jac(i)%off(:,:,k), node(kth_nghbr)%du)

    end do

!   Update du by the GS relaxation (with a relaxation parameter, omega)

     b = matmul( jac(i)%diag_inverse, b ) - node(i)%du
	 linear_res_norm(:,1) = linear_res_norm(:,1) + abs( b(:) )
     node(i)%du = node(i)%du + omega * b

   end do nodes3
   
   linear_res_norm(:,1) = linear_res_norm(:,1) / real(nnodes, p2)

!    End of 1 Sequential Gauss-Seidel Relaxation(sweep)
!----------------------------------------------------

   !Skip the first relaxation
    if (ii==1) then

     linear_res_norm_initial = linear_res_norm
	 write(1000,'(a,i10,a,es12.5)') " relax ", ii, " max(L1 norm) = ", maxval(linear_res_norm(:,1))

   !Check convergence from the second sweep.
    else

      roc = maxval(linear_res_norm(:,1)/linear_res_norm_initial(:,1))

     !write(*,'(a,i10,a,3es12.5,2f10.3,2f10.3)') " relax ", ii, &
     !      " max(L1 norm) = ", linear_res_norm(:,1), roc, omega
	      write(1000,'(a,i10,a,3es12.5,2f10.3,2f10.3)') " relax ", ii, &
           " max(L1 norm) = ", linear_res_norm(:,1), roc, omega

       !Tolerance met. Good!
       if (roc < tolerance_linear .or. &
           maxval(linear_res_norm(:,1)) < 1.0e-17 ) then
       write(1000,*) " Tolerance met. Exit GS relaxation. Total sweeps = ", ii, tolerance_linear
       sweeps_actual = ii
       exit relax

       !Maximum sweep reached... Not converged...
       else

        if (ii == sweeps) then
         write(1000,*) " Tolerance not met... sweeps = ", sweeps
         sweeps_actual = sweeps
        endif

       endif

    endif

!----------------------------------------------------------------------------------------

    linear_res_norm_previous = linear_res_norm

  end do relax

!    End of Linear Relaxation (Sweep)
!--------------------------------------------------------------------------------


  write(1000,*)

 end subroutine gs_sequential2
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------



!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!* This subroutine relaxes the linear system, Jac*du = -Res,
!* by the sequential Symmetric-Gauss-Seidel relaxation scheme.
!*
!* ------------------------------------------------------------------------------
!*  Input:         jac = Jacobian matirx (left hand side)
!*         node(:)%res = Residual (right hand side)
!*              sweeps = The number of GS sweeps to be performed
!*    tolerance_linear = Tolerance on the linear residual
!*
!*
!* Output: node(:)%du  = Correction
!*       sweeps_actual = Number of sweeps performed
!*                 roc = Rate of convergence (the ratio of the final linear residual
!*                       to the initial linear residual), showing how much it is
!*                       reduced. Tolerance is typically 1 order of magnitude.
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine sgs_sequential(sweeps_actual,roc)

 use edu2d_constants   , only : p2, zero
 use edu2d_my_main_data, only : nnodes, node, jac, sweeps, nq, &
                                tolerance_linear, i_iteration

 implicit none

 integer , intent(out)      :: sweeps_actual
 real(p2), intent(out)      :: roc

!Local variables
 real(p2), dimension(nq,nq) :: inverse
 real(p2), dimension(nq)    :: b, x
 integer                    :: i, k, idestat

 integer                    :: ii, kth_nghbr
 real(p2), dimension(nq,3)  :: linear_res_norm !Residual norms(L1,L2,Linf) for the linear system
 real(p2), dimension(nq,3)  :: linear_res_norm_initial
 real(p2), dimension(nq,3)  :: linear_res_norm_previous

 write(1000,*) "DC Iteration = ", i_iteration

 sweeps_actual = 0

!--------------------------------------------------------------------------------
! 1. Initialize the correction

  nodes1 : do i = 1, nnodes
   node(i)%du = zero
  end do nodes1
!--------------------------------------------------------------------------------
! 2. Linear Relaxation (Sweep)

  relax : do ii = 1, sweeps

!----------------------------------------------------
!   Forward Sweep
   linear_res_norm(:,1) = zero

   nodes2 : do i = 1, nnodes

!   Form the right hand side of GS: [ sum( off_diagonal_block*du ) - residual ]

     b = node(i)%res ! Residual has already been given the minus sign.

    do k = 1, node(i)%nnghbrs

     kth_nghbr = node(i)%nghbr(k)
     b = b - matmul(jac(i)%off(:,:,k), node(kth_nghbr)%du)

    end do

!   Update du by the GS relaxation (with a relaxation parameter, omega)

     b = matmul( jac(i)%diag_inverse, b ) - node(i)%du
	 !linear_res_norm(:,1) = linear_res_norm(:,1) + abs( b(:) )
     node(i)%du = node(i)%du + b

   end do nodes2
!
!   Backward sweep
! 
   nodes3 : do i = nnodes, 1, -1

!   Form the right hand side of GS: [ sum( off_diagonal_block*du ) - residual ]

     b = node(i)%res ! Residual has already been given the minus sign.

    do k = 1, node(i)%nnghbrs

     kth_nghbr = node(i)%nghbr(k)
     b = b - matmul(jac(i)%off(:,:,k), node(kth_nghbr)%du)

    end do

!   Update du by the GS relaxation (with a relaxation parameter, omega)

     b = matmul( jac(i)%diag_inverse, b ) - node(i)%du
	 linear_res_norm(:,1) = linear_res_norm(:,1) + abs( b(:) )
     node(i)%du = node(i)%du + b

   end do nodes3
   
   linear_res_norm(:,1) = linear_res_norm(:,1) / real(nnodes, p2)

!    End of Sequential Symmetric-Gauss-Seidel Relaxation(sweep)
!----------------------------------------------------

   !Skip the first relaxation
    if (ii==1) then

     linear_res_norm_initial = linear_res_norm
	 write(1000,'(a,i10,a,es12.5)') " relax ", ii, " max(L1 norm) = ", maxval(linear_res_norm(:,1))

   !Check convergence from the second sweep.
    else

      roc = maxval(linear_res_norm(:,1)/linear_res_norm_initial(:,1))

     !write(*,'(a,i10,a,3es12.5,2f10.3,2f10.3)') " relax ", ii, &
     !      " max(L1 norm) = ", linear_res_norm(:,1), roc, omega
	      write(1000,'(a,i10,a,3es12.5,2f10.3,2f10.3)') " relax ", ii, &
           " max(L1 norm) = ", linear_res_norm(:,1), roc, 1.0

       !Tolerance met. Good!
       if (roc < tolerance_linear .or. &
           maxval(linear_res_norm(:,1)) < 1.0e-17 ) then
       write(1000,*) " Tolerance met. Exit GS relaxation. Total sweeps = ", ii, tolerance_linear
       sweeps_actual = ii
       exit relax

       !Maximum sweep reached... Not converged...
       else

        if (ii == sweeps) then
         write(1000,*) " Tolerance not met... sweeps = ", sweeps
         sweeps_actual = sweeps
        endif

       endif

    endif

!----------------------------------------------------------------------------------------

    linear_res_norm_previous = linear_res_norm

  end do relax

!    End of Linear Relaxation (Sweep)
!--------------------------------------------------------------------------------


  write(1000,*)

 end subroutine sgs_sequential
!--------------------------------------------------------------------------------
 end module edu2d_euler_linear_solve
