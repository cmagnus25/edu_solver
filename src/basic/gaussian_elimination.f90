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
 module gaussian_elimination
 
 private

 public :: gewp_solve

 contains

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
 end module gaussian_elimination