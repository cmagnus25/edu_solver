module vector_operations

 private

 public :: elltwo_norm
 public :: rms
 public :: length

 contains

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