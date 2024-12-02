module vector_operations

 private

 public :: w2u
 public :: u2w
 public :: elltwo_norm
 public :: rms
 public :: length
 public :: vector_res
 public :: vector_u
 public :: vector2du
 public :: vector_dt_term

 contains
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
 use input_parameter   , only : gamma

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
 use input_parameter   , only : gamma

 implicit none

 real(p2), dimension(4), intent(in) :: u

!Local variables
 real(p2), dimension(4)             :: w !output

  w(1) = u(1)
  w(2) = u(2)/u(1)
  w(3) = u(3)/u(1)
  w(4) = (gamma-one)*( u(4) - half*w(1)*(w(2)*w(2)+w(3)*w(3)) )

 end function u2w
!********************************************************************************
 function vector_res(nnodes,nq) result(vector)

 use edu2d_constants   , only : p2
 use edu2d_my_main_data, only : node

 implicit none

!Input
 integer          , intent( in) :: nnodes, nq

!Output
 real(p2), dimension(nnodes*nq) :: vector

!local variables
 integer                        :: i, k

  do i = 1, nnodes
   do k = 1, nq
    vector( nq*(i-1) + k ) = node(i)%res(k)
   end do
  end do  

 end function vector_res
!********************************************************************************

!********************************************************************************
!* This function returns a global vector of residual.
!********************************************************************************
 subroutine vector_u(vector)

 use edu2d_constants   , only : p2
 use edu2d_my_main_data, only : node, nnodes, nq

 implicit none

!Output
 real(p2), dimension(nnodes*nq), intent(out) :: vector

!local variables
 integer                                     :: i, k

  do i = 1, nnodes
   do k = 1, nq
    vector( nq*(i-1) + k ) = node(i)%u(k)
   end do
  end do  

 end subroutine vector_u
!********************************************************************************


!********************************************************************************
!* This function returns a global vector of residual.
!********************************************************************************
 subroutine vector2du(vector,nnodes,nq)

 use edu2d_constants   , only : p2
 use edu2d_my_main_data, only : node

 implicit none

!Input
 integer                       , intent( in) :: nnodes, nq
 real(p2), dimension(nnodes*nq), intent( in) :: vector

!local variables
 integer                                     :: i, k

  do i = 1, nnodes
   do k = 1, nq
    node(i)%du(k) = vector( nq*(i-1) + k )
   end do
  end do  

 end subroutine vector2du
 !********************************************************************************
 subroutine vector_dt_term(vector)

 use edu2d_constants   , only : p2
 use edu2d_my_main_data, only : node, nnodes, nq, CFL_frechet, CFL

 implicit none

!Output
 real(p2), dimension(nnodes*nq), intent(out) :: vector

!local variables
 integer                                     :: i, k
 real(p2) :: temp

  CFL_frechet = CFL

  do i = 1, nnodes
   temp = node(i)%vol / ( CFL_frechet*node(i)%dt)

   do k = 1, nq
    vector( nq*(i-1) + k ) = temp
   end do
  end do  

 end subroutine vector_dt_term
!********************************************************************************
 function elltwo_norm(vec,n)

 use edu2d_constants, only : p2, zero

 implicit none

!Input
 integer               , intent(in) :: n
 real(p2), dimension(n), intent(in) :: vec

!Output
 real(p2)                           :: elltwo_norm

!Local variable
 integer :: i

   elltwo_norm = zero

  do i = 1, n
   elltwo_norm = elltwo_norm + abs(vec(i))**2
  end do

   elltwo_norm = sqrt( elltwo_norm / real(n,p2) )

 end function elltwo_norm
!********************************************************************************
 function rms(vec,n)

 use edu2d_constants, only : p2, zero

 implicit none

!Input
 integer               , intent(in) :: n
 real(p2), dimension(n), intent(in) :: vec
!Output
 real(p2)                           :: rms
!Local variable
 integer :: i

   rms = zero

  do i = 1, n
   rms = rms + vec(i)**2
  end do

   rms = sqrt( rms / real(n,p2) )

 end function rms
!********************************************************************************
 function length(vec,n)

 use edu2d_constants, only : p2, zero

 implicit none

!Input
 integer               , intent(in) :: n
 real(p2), dimension(n), intent(in) :: vec

!Output
 real(p2)                           :: length

!Local variable
 integer :: i

   length = zero

  do i = 1, n
   length = length + vec(i)**2
  end do

   length = sqrt( length )

 end function length
!********************************************************************************
end module vector_operations