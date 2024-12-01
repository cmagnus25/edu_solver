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
 end module edu2d_euler_linear_solve
