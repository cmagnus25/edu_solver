!********************************************************************************
!* This subroutine performs preconditioning by the implicit DC solver.
!*
!*
!* linear_solve with right hand side "r"
!*
!* Written by Katate Masatsuka (http://www.cfdbooks.com)
!********************************************************************************
module precond

 private
 
 public :: preconditioner
   
 contains
 
 subroutine preconditioner(r,ndim, p,istat)

 use edu2d_constants   , only : p2
 use edu2d_my_main_data, only : nq, nnodes, node, sweeps_actual_gcr
 use edu2d_euler_linear_solve, only : smooth

 implicit none

!Input
 integer                  , intent( in) :: ndim
 real(p2), dimension(ndim), intent( in) :: r
!Output
 real(p2), dimension(ndim), intent(out) :: p
 integer ,                  intent(out) :: istat

!Local variables
 integer                                :: i, k, sweeps_actual
 real(p2)                               :: roc

!--------------------------------------------------------------------------------
! Save the current residual

   do i = 1, nnodes
     node(i)%r_temp = node(i)%res
   end do

!--------------------------------------------------------------------------------
! Store the input residual into the residual vector.

  do i = 1, nnodes
   do k = 1, nq
    node(i)%res(k) = r( nq*(i-1) + k )
   end do
  end do

!--------------------------------------------------------------------------------
! Linear relaxation: Relax the linear system of Jac*dU = r,
!                    which is the implciit Defect Correction solver.

! Relax Jac*du = r by sequential/multicolor GS.
  call smooth(sweeps_actual,roc) !This will compute node(:)%du(:).
  write(2000,*) "preconditioning:", sweeps_actual,":",roc
  write(*,'(27x,a32,i5,a,es8.2,a5,es10.2,a5,es10.2)') " - preconditioner: Smooth(sweeps:cr)=", &
  sweeps_actual,":",roc
  sweeps_actual_gcr = sweeps_actual_gcr + sweeps_actual

!--------------------------------------------------------------------------------
! Store the solution into p(:):
! ->   p(:) = -(Jac)^{-1}*r, where the inverse is approximate.

  do i = 1, nnodes
   do k = 1, nq
    p( nq*(i-1) + k ) = node(i)%du(k)
   end do
  end do

!--------------------------------------------------------------------------------
! Copy back the current residual.

   do i = 1, nnodes
     node(i)%res = node(i)%r_temp
   end do

!--------------------------------------------------------------------------------
  if (roc > 1.0_p2) then
   istat = 1
  else
   istat = 0
  endif

 end subroutine preconditioner
!********************************************************************************
end module precond