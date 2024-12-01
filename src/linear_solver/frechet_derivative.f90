!********************************************************************************
!* This subroutine computes the Frechet derivative
!*
!* Ap = dRes/dU*p ~ [Res(u+epsilon*p) - Res(u)]/epsilon
!* Note: compute_residual_ncfv computes %res = -Res.
!*       So, don't forget the minus sign to compute Ap.
!*
!* Written by Katate Masatsuka (http://www.cfdbooks.com)
!********************************************************************************
module frechet_derivative

private

public :: compute_Ap

contains

 function compute_Ap(p,r0,ndim) result(Ap)

 use edu2d_constants               , only : p2, one, zero
 use edu2d_my_main_data            , only : nq, nnodes, node, my_eps
 use residual , only : compute_residual_ncfv
 use vector_operations , only : length, vector_u, rms, vector_dt_term, u2w

 implicit none

!Input
 integer                  , intent( in) :: ndim
 real(p2), dimension(ndim), intent( in) :: p, r0
!Output
 real(p2), dimension(ndim) :: Ap

!Local variables
 real(p2), dimension(ndim) :: r_ep, v_over_dt
 real(p2)                  :: eps, length_p
 integer                   :: i, k

 Ap = zero

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Frechet derivative: Newton's Method

!--------------------------------------------------------------------------------
! Save the current solution and residual
! Need to store temporarily primitive variables

   do i = 1, nnodes
     node(i)%u_temp = node(i)%u
	 node(i)%w_temp = node(i)%w
     node(i)%r_temp = node(i)%res ! This is -Res.
   end do

!--------------------------------------------------------------------------------
! Store (u + epsilon*p) into the solution vector.

   length_p = length(p,ndim)

   call vector_u(r_ep) !Store the solution temporarily in r_ep(:).
   eps = max( one, rms(r_ep,ndim) )*sqrt(my_eps)

   do i = 1, nnodes
    do k = 1, nq
     node(i)%u(k) = node(i)%u(k) + eps* p( nq*(i-1) + k )/max(my_eps,length_p)
    end do
   end do
!--------------------------------------------------------------------------------
!  Need to recompute primitive variables 
   do i = 1, nnodes
     node(i)%w = u2w(node(i)%u)
   end do
!--------------------------------------------------------------------------------
! Compute the residual vector with u + epsilon*p.

   call compute_residual_ncfv

!--------------------------------------------------------------------------------
! Store the computed residual into the residual vector, r_ep.

   do i = 1, nnodes
    do k = 1, nq
     r_ep( nq*(i-1) + k ) = node(i)%res(k) ! This is -Res.
    end do
   end do

!--------------------------------------------------------------------------------
! Compute the Frechet derivative (finite-difference approximation).

   Ap = ( (-r_ep) - (-r0) ) / eps * length_p ! = [Res(u+epsilon*p) - Res(u)]/epsilon

!--------------------------------------------------------------------------------
! Compute the pseudo-time term

   call vector_dt_term(v_over_dt)

   do i = 1, ndim
    Ap(i) = v_over_dt(i)*p(i) + Ap(i)  ! Ap <- (V/dt + A)*p
   end do

!--------------------------------------------------------------------------------
! Copy back the current solution and residual

   do i = 1, nnodes
     node(i)%u   = node(i)%u_temp
	 node(i)%w   = node(i)%w_temp
     node(i)%res = node(i)%r_temp ! This is -Res.
   end do

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

 end function compute_Ap
!********************************************************************************
end module frechet_derivative