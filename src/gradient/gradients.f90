!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!* Below, you find the following subroutines associated with LSQ gradients:
!*
!*  - compute_lsq_coeff_nc : Compute and store LSQ gradient coefficients at node
!*  - check_lsq_coeff_nc   : Check the computed LSQ coefficients
!*  - compute_gradient_nc  : Compute gradients at nodes
!*
!*  - lsq_gradients_nc     : Compute gradients (loop over neighbros) at a node
!*  - lsq_gradients2_nc    : Compute gradients (including neighbors of neighbors) at a node
!*  - lsq01_2x2_coeff_nc   : Compute linear LSQ coefficients at a node
!*  - lsq02_5x5_coeff2_nc  : Compute quadratic LSQ coefficients at all nodes
!*  - lsq_weight           : Compute LSQ weights
!
!********************************************************************************
!********************************************************************************

!********************************************************************************
!* This subroutine computes the coefficients for linear and quadratic LSQ
!* gradeints. 
!*
!* Note: Quadratic LSQ method is implemented in two steps as described in
!*       Nishikawa, JCP2014v273pp287-309 for details, which is available at 
!* http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2014v273pp287-309_preprint.pdf.
!*       This two-step method is useful in a paralell environment: it takes into
!*       account neighbors of neighbors without accessing neighbors of neighbors.
!*       So, it can be implemented only with edge-connected-neighbor information.
!*       In other words, each step is compact.
!*
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* This is Version 0 (July 2015).
!* This F90 code is written and made available for an educational purpose.
!* This file may be updated in future.
!*
!********************************************************************************

 module gradients_lsq

 private

 public :: compute_lsq_coeff_nc
 public ::   check_lsq_coeff_nc
 public :: compute_gradient_nc

 contains

 subroutine compute_lsq_coeff_nc

 use edu2d_my_main_data , only : nnodes, node, nq
 use edu2d_my_allocation, only : my_alloc_p2_ptr, my_alloc_p2_matrix_ptr

 integer :: i, in, ell, ii, k

  write(*,*)
  write(*,*) "Constructing LSQ coefficients..."

! 1. Coefficients for the linear LSQ gradients

  write(*,*) "---(1) Constructing Linear LSQ coefficients..."

  nodes : do i = 1, nnodes

     call my_alloc_p2_ptr(node(i)%lsq2x2_cx,node(i)%nnghbrs)
     call my_alloc_p2_ptr(node(i)%lsq2x2_cy,node(i)%nnghbrs)
     call lsq01_2x2_coeff_nc(i)

  end do nodes

! 2. Coefficients for the quadratic LSQ gradients (two-step method)

  write(*,*) "---(2) Constructing Quadratic LSQ coefficients..."

  do i = 1, nnodes

        ii = 0 
     nghbr : do k = 1, node(i)%nnghbrs
        in = node(i)%nghbr(k)
      nghbr_nghbr : do ell = 1, node(in)%nnghbrs
        ii = ii + 1
      end do nghbr_nghbr
     end do nghbr

     call my_alloc_p2_ptr(node(i)%lsq5x5_cx, ii)
     call my_alloc_p2_ptr(node(i)%lsq5x5_cy, ii)
     call my_alloc_p2_ptr(node(i)%dx,node(i)%nnghbrs)
     call my_alloc_p2_ptr(node(i)%dy,node(i)%nnghbrs)
     call my_alloc_p2_matrix_ptr(node(i)%dw, nq,node(i)%nnghbrs)

  end do

     call lsq02_5x5_coeff2_nc

 end subroutine compute_lsq_coeff_nc

!********************************************************************************
!* This subroutine verifies the implementation of LSQ gradients.
!*
!* 1. Check if the linear LSQ gradients are exact for linear functions.
!* 2. Check if the quadratic LSQ gradients are exact for quadratic functions.
!*
!* Note: Here, we use only the first component of u=(u1,u2,u3), i.e., ivar=1.
!*
!********************************************************************************
 subroutine check_lsq_coeff_nc

 use edu2d_constants   , only : p2, one, two
 use edu2d_my_main_data, only : nnodes, node

 integer       :: i, ix, iy, ivar
 character(80) :: grad_type_temp
 real(p2)      :: error_max_wx, error_max_wy, x, y
 real(p2)      :: x_max_wx, y_max_wx, x_max_wy, y_max_wy, wx, wxe, wy, wye
 real(p2)      :: a0, a1, a2, a3, a4, a5

  ix = 1
  iy = 2

! We only use w(1) for this test.
  ivar = 1
 
!---------------------------------------------------------------------
! 1. Check linear LSQ gradients
!---------------------------------------------------------------------
  write(*,*)
  write(*,*) "---------------------------------------------------------"
  write(*,*) "---------------------------------------------------------"
  write(*,*) "- Checking Linear LSQ gradients..."

!  (1). Store a linear function in w(ivar) = x + 2*y.
!       So the exact gradient is grad(w(ivar)) = (1,2).

   write(*,*) "- Storing a linear function values..."
   do i = 1, nnodes
    x = node(i)%x
    y = node(i)%y
    node(i)%w(ivar) = one*x + two*y
   end do

!  (2). Compute the gradient by linear LSQ

   write(*,*) "- Computing linear LSQ gradients.."
   grad_type_temp = 'linear'
   call compute_gradient_nc(ivar,grad_type_temp)

!  (3). Compute the relative errors (L_infinity)

   write(*,*) "- Computing the relative errors (L_infinity).."
   error_max_wx = -one
   error_max_wy = -one
   do i = 1, nnodes
    error_max_wx = max( abs( node(i)%gradw(ivar,ix) - one )/one, error_max_wx )
    error_max_wy = max( abs( node(i)%gradw(ivar,iy) - two )/two, error_max_wy )
   end do

  write(*,*) " Max relative error in wx = ", error_max_wx
  write(*,*) " Max relative error in wy = ", error_max_wy
  write(*,*) "---------------------------------------------------------"
  write(*,*) "---------------------------------------------------------"


!---------------------------------------------------------------------
! 2. Check quadratic LSQ gradients
!---------------------------------------------------------------------
  write(*,*)
  write(*,*) "---------------------------------------------------------"
  write(*,*) "---------------------------------------------------------"
  write(*,*) "- Checking Quadratic LSQ gradients..."

!  (1). Store a quadratic function in w(ivar) = a0 + a1*x + a2*y + a3*x**2 + a4*x*y + a5*y**2
!       So the exact gradient is grad(w(ivar)) = (a1+2*a3*x+a4*y, a2+2*a5*y+a4*x)

   a0 =     21.21_p2
   a1 =     1.000_p2
   a2 = -   1.970_p2
   a3 =   280.400_p2
   a4 = -2129.710_p2
   a5 =   170.999_p2

   write(*,*) "- Storing a quadratic function values..."
   do i = 1, nnodes
    x = node(i)%x
    y = node(i)%y
    node(i)%w(ivar) = a0 + a1*x + a2*y + a3*x**2 + a4*x*y + a5*y**2
   end do

!  (2). Compute the gradient by linear LSQ

   write(*,*) "- Computing quadratic LSQ gradients.."
   grad_type_temp = 'quadratic2'
   call compute_gradient_nc(ivar,grad_type_temp)

!  (3). Compute the relative errors (L_infinity)

   write(*,*) "- Computing the relative errors (L_infinity).."
   error_max_wx = -one
   error_max_wy = -one
   do i = 1, nnodes
    x = node(i)%x
    y = node(i)%y

    if ( abs( node(i)%gradw(ivar,ix) - (a1+2.0_p2*a3*x+a4*y) )/(a1+2.0_p2*a3*x+a4*y) >  error_max_wx ) then
      wx  = node(i)%gradw(ivar,ix)
      wxe = a1+2.0_p2*a3*x+a4*y
      error_max_wx = abs( wx - wxe )/wxe
      x_max_wx = x
      y_max_wx = y
    endif

    if ( abs( node(i)%gradw(ivar,iy) - (a2+2.0_p2*a5*y+a4*x) )/(a2+2.0_p2*a5*y+a4*x) >  error_max_wy ) then
      wy  = node(i)%gradw(ivar,iy)
      wye = a2+2.0_p2*a5*y+a4*x
      error_max_wy = abs( wy - wye )/wye
      x_max_wy = x
      y_max_wy = y
    endif

   end do

  write(*,'(a,es20.3,a,2es12.5)') " Max relative error in wx = ", error_max_wx, " at (x,y) = ", x_max_wx, y_max_wx
  write(*,'(a,es20.10,a,es20.10)')  "   At this location, LSQ ux = ", wx, ": Exact ux = ", wxe
  write(*,'(a,es20.3,a,2es12.5)') " Max relative error in wy = ", error_max_wy, " at (x,y) = ", x_max_wy, y_max_wy
  write(*,'(a,es20.10,a,es20.10)')  "   At this location, LSQ uy = ", wy, ": Exact uy = ", wye
  write(*,*) "---------------------------------------------------------"
  write(*,*) "---------------------------------------------------------"
  write(*,*)


 end subroutine check_lsq_coeff_nc

!********************************************************************************
!* This subroutine computes gradients at nodes for the variable u(ivar),
!* where ivar = 1,2,3, ..., or nq.
!*
!* ------------------------------------------------------------------------------
!*  Input: node(:)%u(ivar)
!*
!* Output: node(i)%gradu(ivar,1:2) = ( du(ivar)/dx, du(ivar)/dy )
!* ------------------------------------------------------------------------------
!********************************************************************************
 subroutine compute_gradient_nc(ivar,grad_type)

 use edu2d_my_main_data, only : node, nnodes

 integer, intent(in) :: ivar

 integer       :: i, k, in
 character(80) :: grad_type

  if (trim(grad_type) == "none") return

!-------------------------------------------------
!  Two-step quadratic LSQ 5x5 system
!  Note: See Nishikawa, JCP2014v273pp287-309 for details, which is available at 
!        http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2014v273pp287-309_preprint.pdf.

   if (trim(grad_type) == "quadratic2") then

!  Perform Step 1 as below (before actually compute the gradient).

      do i = 1, nnodes

       nghbr0 : do k = 1, node(i)%nnghbrs
                   in      = node(i )%nghbr(k)
        node(i)%dx(k)      = node(in)%x       - node(i)%x
        node(i)%dy(k)      = node(in)%y       - node(i)%y
        node(i)%dw(ivar,k) = node(in)%w(ivar) - node(i)%w(ivar)
       end do nghbr0
      end do

   endif
!-------------------------------------------------

!------------------------------------------------------------
!------------------------------------------------------------
!-- Compute LSQ Gradients at all nodes.
!------------------------------------------------------------
!------------------------------------------------------------

  nodes : do i = 1, nnodes

  !-------------------------------------------------
  ! Linear LSQ 2x2 system
    if (trim(grad_type) == "linear") then

       call lsq_gradients_nc(i,ivar)

  !-------------------------------------------------
  ! Two-step quadratic LSQ 5x5 system
  !  Note: See Nishikawa, JCP2014v273pp287-309 for details, which is available at 
  !        http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2014v273pp287-309_preprint.pdf.
    elseif (trim(grad_type) == "quadratic2") then

      call lsq_gradients2_nc(i,ivar)

  !-------------------------------------------------
    else

     write(*,*) " Invalid input value -> ", trim(grad_type)
     stop

    endif
  !-------------------------------------------------

   end do nodes

 end subroutine compute_gradient_nc


!********************************************************************************
!* Compute the gradient, (wx,wy), for the variable u by Linear LSQ.
!*
!* ------------------------------------------------------------------------------
!*  Input:            inode = Node number at which the gradient is computed.
!*                     ivar =   Variable for which the gradient is computed.
!*          node(:)%w(ivar) = Solution at nearby nodes.
!*
!* Output:  node(inode)%gradu = gradient of the requested variable
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine lsq_gradients_nc(inode,ivar)

 use edu2d_my_main_data           , only : node
 use edu2d_constants              , only : p2, zero

 implicit none

 integer, intent(in) :: inode, ivar

!Local variables
 integer  :: in, inghbr
 integer  :: ix, iy
 real(p2) :: da, ax, ay

  ix = 1 
  iy = 2
  ax = zero
  ay = zero

!   Loop over neighbors

     do in = 1, node(inode)%nnghbrs
       inghbr = node(inode)%nghbr(in)

          da = node(inghbr)%w(ivar) - node(inode)%w(ivar)
 
      ax = ax + node(inode)%lsq2x2_cx(in)*da
      ay = ay + node(inode)%lsq2x2_cy(in)*da

     end do

      node(inode)%gradw(ivar,ix) = ax  !<-- dw(ivar)/dx
      node(inode)%gradw(ivar,iy) = ay  !<-- dw(ivar)/dy

 end subroutine lsq_gradients_nc
!--------------------------------------------------------------------------------


!********************************************************************************
!* Compute the gradient, (wx,wy), for the variable u by Quadratic LSQ.
!*
!* ------------------------------------------------------------------------------
!*  Input:            inode = Node number at which the gradient is computed.
!*                     ivar =   Variable for which the gradient is computed.
!*          node(:)%w(ivar) = Solution at nearby nodes.
!*
!* Output:  node(inode)%gradu = gradient of the requested variable
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine lsq_gradients2_nc(inode,ivar)

 use edu2d_my_main_data, only : node
 use edu2d_constants   , only : p2, zero

 implicit none

 integer, intent(in) :: inode, ivar

!Local variables
 integer  :: in
 integer  :: ix, iy, ii, ell, k
 real(p2) :: da, ax, ay

   ix = 1 
   iy = 2
   ax = zero
   ay = zero

!   Loop over neighbors

       ii = 0

     nghbr : do k = 1, node(inode)%nnghbrs
              in = node(inode)%nghbr(k)

      nghbr_nghbr : do ell = 1, node(in)%nnghbrs

       da = node(in)%w(ivar) - node(inode)%w(ivar) + node(in)%dw(ivar,ell)

       if ( node(in)%nghbr(ell) == inode ) then
        da = node(in)%w(ivar) - node(inode)%w(ivar)
       endif

       ii = ii + 1

       ax = ax + node(inode)%lsq5x5_cx(ii)*da
       ay = ay + node(inode)%lsq5x5_cy(ii)*da

      end do nghbr_nghbr

     end do nghbr

      node(inode)%gradw(ivar,ix) = ax  !<-- dw(ivar)/dx
      node(inode)%gradw(ivar,iy) = ay  !<-- dw(ivar)/dy

 end subroutine lsq_gradients2_nc
!--------------------------------------------------------------------------------



!********************************************************************************
!* --- LSQ Coefficients for 2x2 Linear Least-Squares Gradient Reconstruction ---
!*
!* ------------------------------------------------------------------------------
!*  Input:  inode = node number at which the gradient is computed.
!*
!* Output:  node(inode)%lsq2x2_cx(:)
!*          node(inode)%lsq2x2_cx(:)
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine lsq01_2x2_coeff_nc(inode)

 use edu2d_my_main_data           , only : node
 use edu2d_constants              , only : p2, zero

 implicit none

 integer, intent(in) :: inode
!Local variables
 real(p2) :: a(2,2), dx, dy, det, w2, w2dvar
 integer  :: k, inghbr, ix=1,iy=2
 real(p2), dimension(2,2) :: local_lsq_inverse

   a = zero

!  Loop over the neighbor nodes.
   do k = 1, node(inode)%nnghbrs
    inghbr = node(inode)%nghbr(k)

      dx = node(inghbr)%x - node(inode)%x
      dy = node(inghbr)%y - node(inode)%y

      w2 = lsq_weight(dx, dy)**2

      a(1,1) = a(1,1) + w2 * dx*dx
      a(1,2) = a(1,2) + w2 * dx*dy

      a(2,1) = a(2,1) + w2 * dx*dy
      a(2,2) = a(2,2) + w2 * dy*dy

   end do

    det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
    if (abs(det) < 1.0e-14) write(*,*) " Singular: LSQ det = ", det, " i=",inode

! OK, invert and store the inverse matrix:

     local_lsq_inverse(1,1) =  a(2,2)/det
     local_lsq_inverse(1,2) = -a(2,1)/det
     local_lsq_inverse(2,1) = -a(1,2)/det
     local_lsq_inverse(2,2) =  a(1,1)/det

!  Now compute the coefficients for neighbors.

     nghbr : do k = 1, node(inode)%nnghbrs
       inghbr = node(inode)%nghbr(k)

        dx = node(inghbr)%x - node(inode)%x
        dy = node(inghbr)%y - node(inode)%y

      w2dvar = lsq_weight(dx, dy)**2

      node(inode)%lsq2x2_cx(k)  = local_lsq_inverse(ix,1)*w2dvar*dx &
                                + local_lsq_inverse(ix,2)*w2dvar*dy

      node(inode)%lsq2x2_cy(k)  = local_lsq_inverse(iy,1)*w2dvar*dx &
                                + local_lsq_inverse(iy,2)*w2dvar*dy

     end do nghbr

 end subroutine lsq01_2x2_coeff_nc
!********************************************************************************
!*
!********************************************************************************


!********************************************************************************
!* --- LSQ Coefficients for 5x5 Quadratic Least-Squares Gradient Reconstruction ---
!*
!* Note: See Nishikawa, JCP2014v273pp287-309 for details, which is available at 
!*
!* http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2014v273pp287-309_preprint.pdf.
!*
!* ------------------------------------------------------------------------------
!*  Input:
!*
!* Output:  node(:)%lsq5x5_cx(ii)
!*          node(:)%lsq5x5_cx(ii)
!*
!* Note: This subroutine computes the LSQ coefficeints at all nodes.
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine lsq02_5x5_coeff2_nc

 use edu2d_my_main_data           , only : node, nnodes
 use edu2d_constants              , only : p2, zero, half
 use gaussian_elimination         , only : gewp_solve

 implicit none

!Local variables
 real(p2) :: a(5,5), ainv(5,5), dx, dy
 real(p2) :: dummy1(5), dummy2(5)
 real(p2) :: w2
 integer  :: istat, ix=1, iy=2

 integer :: i, k, ell, in, ii

! Step 1

  node1 : do i = 1, nnodes
   nghbr : do k = 1, node(i)%nnghbrs

               in   = node(i)%nghbr(k)
    node(i)%dx(k)   = node(in)%x - node(i)%x
    node(i)%dy(k)   = node(in)%y - node(i)%y

   end do nghbr
  end do node1

! Step 2

  node2 : do i = 1, nnodes

     a = zero

!  Get dx, dy, and dw

   nghbr2 : do k = 1, node(i)%nnghbrs

              in = node(i)%nghbr(k)

    nghbr_nghbr : do ell = 1, node(in)%nnghbrs

     dx = node(in)%x - node(i)%x + node(in)%dx(ell)
     dy = node(in)%y - node(i)%y + node(in)%dy(ell)

     if ( node(in)%nghbr(ell) == i ) then

      dx = node(in)%x - node(i)%x
      dy = node(in)%y - node(i)%y

      if ( abs(dx) + abs(dy) < 1.0e-13_p2 ) then
       write(*,*) " Zero distance found at lsq02_5x5_coeff2_nc..."
       write(*,*) "    dx = ", dx
       write(*,*) "    dy = ", dy
       write(*,*) "- Centered node = ", i
       write(*,*) "          (x,y) = ", node(in)%x, node(in)%y
       write(*,*) "- Neighbor node = ", in
       write(*,*) "          (x,y) = ", node(in)%x, node(in)%y
       stop
      endif

     endif


        w2 = lsq_weight(dx, dy)**2

      a(1,1) = a(1,1) + w2 * dx         *dx
      a(1,2) = a(1,2) + w2 * dx         *dy
      a(1,3) = a(1,3) + w2 * dx         *dx*dx * half
      a(1,4) = a(1,4) + w2 * dx         *dx*dy
      a(1,5) = a(1,5) + w2 * dx         *dy*dy * half

!     a(2,1) = a(2,1) + w2 * dy         *dx
      a(2,2) = a(2,2) + w2 * dy         *dy
      a(2,3) = a(2,3) + w2 * dy         *dx*dx * half
      a(2,4) = a(2,4) + w2 * dy         *dx*dy
      a(2,5) = a(2,5) + w2 * dy         *dy*dy * half

!     a(3,1) = a(3,1) + w2 * half*dx*dx *dx
!     a(3,2) = a(3,2) + w2 * half*dx*dx *dy
      a(3,3) = a(3,3) + w2 * half*dx*dx *dx*dx * half
      a(3,4) = a(3,4) + w2 * half*dx*dx *dx*dy
      a(3,5) = a(3,5) + w2 * half*dx*dx *dy*dy * half

!     a(4,1) = a(4,1) + w2 *      dx*dy *dx
!     a(4,2) = a(4,2) + w2 *      dx*dy *dy
!     a(4,3) = a(4,3) + w2 *      dx*dy *dx*dx * half
      a(4,4) = a(4,4) + w2 *      dx*dy *dx*dy
      a(4,5) = a(4,5) + w2 *      dx*dy *dy*dy * half

!     a(5,1) = a(5,1) + w2 * half*dy*dy *dx
!     a(5,2) = a(5,2) + w2 * half*dy*dy *dy
!     a(5,3) = a(5,3) + w2 * half*dy*dy *dx*dx * half
!     a(5,4) = a(5,4) + w2 * half*dy*dy *dx*dy
      a(5,5) = a(5,5) + w2 * half*dy*dy *dy*dy * half

    end do nghbr_nghbr

   end do nghbr2

!   Fill symmetric part

      a(2,1) = a(1,2);
      a(3,1) = a(1,3);  a(3,2) = a(2,3);
      a(4,1) = a(1,4);  a(4,2) = a(2,4); a(4,3) = a(3,4);
      a(5,1) = a(1,5);  a(5,2) = a(2,5); a(5,3) = a(3,5); a(5,4) = a(4,5);

!   Invert the matrix

     dummy1 = zero
     dummy2 = zero
     call gewp_solve(a,dummy1,dummy2,ainv,istat, 5);

     if (istat/=0) write(*,*) "Problem in solving the linear system!: Quadratic_LSJ_Matrix"

!  Now compute the coefficients for neighbors.

      ii = 0

     nghbr3 : do k = 1, node(i)%nnghbrs
                in = node(i)%nghbr(k)

      nghbr_nghbr2 : do ell = 1, node(in)%nnghbrs

       dx = node(in)%x - node(i)%x + node(in)%dx(ell)
       dy = node(in)%y - node(i)%y + node(in)%dy(ell)

       if ( node(in)%nghbr(ell) == i ) then
        dx = node(in)%x - node(i)%x
        dy = node(in)%y - node(i)%y
       endif

       ii = ii + 1
       
       w2 = lsq_weight(dx, dy)**2

 !  Multiply the inverse LSQ matrix to get the coefficients: cx(:) and cy(:):

      node(i)%lsq5x5_cx(ii)  =  ainv(ix,1)*w2*dx             &
                              + ainv(ix,2)*w2*dy             &
                              + ainv(ix,3)*w2*dx*dx * half   &
                              + ainv(ix,4)*w2*dx*dy          &
                              + ainv(ix,5)*w2*dy*dy * half

      node(i)%lsq5x5_cy(ii)  =  ainv(iy,1)*w2*dx             &
                              + ainv(iy,2)*w2*dy             &
                              + ainv(iy,3)*w2*dx*dx * half   &
                              + ainv(iy,4)*w2*dx*dy          &
                              + ainv(iy,5)*w2*dy*dy * half
      end do nghbr_nghbr2

     end do nghbr3

  end do node2

 end subroutine lsq02_5x5_coeff2_nc
!********************************************************************************
!*
!********************************************************************************


!****************************************************************************
!* Compute the LSQ weight
!*
!* Note: The weight computed here is the square of the actual LSQ weight.
!*****************************************************************************
 function lsq_weight(dx, dy)

 use edu2d_constants   , only : p2, one
 use edu2d_my_main_data, only : gradient_weight, gradient_weight_p

 implicit none

!Input
 real(p2), intent(in) :: dx, dy
!Output
 real(p2)             :: lsq_weight
!Local
 real(p2)             :: distance

  if     (trim(gradient_weight) == "none"            ) then

        lsq_weight = one

  elseif (trim(gradient_weight) == "inverse_distance") then

          distance = sqrt(dx*dx + dy*dy)
        lsq_weight = one / distance**gradient_weight_p

  endif

 end function lsq_weight
 !--------------------------------------------------------------------------------
 end module gradients_lsq