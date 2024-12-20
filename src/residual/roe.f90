!--------------------------------------------------------------------------------


!********************************************************************************
!* -- Roe's Flux Function with entropy fix---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!*
!* NOTE: 3D version of this subroutine is available for download at
!*       http://cfdbooks.com/cfdcodes.html
!*
!* ------------------------------------------------------------------------------
!*  Input:   primL(1:4) =  left state (rhoL, uL, vL, pL)
!*           primR(1:4) = right state (rhoR, uR, vR, pR)
!*               njk(2) = Face normal (L -> R). Must be a unit vector.
!*
!* Output:    flux(1:4) = numerical flux
!*                  wsn = half the max wave speed
!*                        (to be used for time step calculations)
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
module roe_scheme
 
 private
 
 public :: roe
 
 contains

 subroutine roe(primL, primR, njk,  flux, wsn)

 use edu2d_constants   , only : p2, zero, one, two, half, fifth
 use input_parameter, only : gamma

 implicit none

!Input:
 real(p2), intent( in) :: primL(4), primR(4) ! [rho, u, v, p]_{L,R}
 real(p2), intent( in) :: njk(2)             ! Face normal, njk=[nx, ny]

!Output
 real(p2), intent(out) :: flux(4)
 real(p2), intent(out) :: wsn

!Local variables
 real(p2) :: nx, ny                  ! Normal vector
 real(p2) :: mx, my                  ! Tangent vector: mx*nx+my*ny = 0
 real(p2) :: uL, uR, vL, vR          ! Velocity components.
 real(p2) :: rhoL, rhoR, pL, pR      ! Primitive variables.
 real(p2) :: unL, unR, umL, umR      ! Normal and tangent velocities
 real(p2) :: aL, aR, HL, HR          ! Speeds of sound.
 real(p2) :: RT,rho,u,v,H,a,un, um   ! Roe-averages
 real(p2) :: drho,dun,dum,dp,LdU(4)  ! Wave strenghs
 real(p2) :: ws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
 real(p2) :: fL(4), fR(4), diss(4)   ! Fluxes ad dissipation term
 real(p2) :: dws(4)                  ! User-specified width for entropy fix
 integer :: i, j

  nx = njk(1)
  ny = njk(2)

!Tangent vector (Do you like it? Actually, Roe flux can be implemented 
! without any tangent vector. See "I do like CFD, VOL.1" for details.)
  mx = -ny
  my =  nx

!Primitive and other variables.
!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
     unL = uL*nx+vL*ny
     umL = uL*mx+vL*my
      pL = primL(4)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL)
!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
     unR = uR*nx+vR*ny
     umR = uR*mx+vR*my
      pR = primR(4)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR)

!First compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
     u = (uL+RT*uR)/(one+RT)
     v = (vL+RT*vR)/(one+RT)
     H = (HL+RT* HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*(u*u+v*v)) )
    un = u*nx+v*ny
    um = u*mx+v*my

!Wave Strengths
   drho = rhoR - rhoL 
     dp =   pR - pL
    dun =  unR - unL
    dum =  umR - umL

  LdU(1) = (dp - rho*a*dun )/(two*a*a)
  LdU(2) = rho*dum
  LdU(3) =  drho - dp/(a*a)
  LdU(4) = (dp + rho*a*dun )/(two*a*a)

!Wave Speed
  ws(1) = abs(un-a)
  ws(2) = abs(un)
  ws(3) = abs(un)
  ws(4) = abs(un+a)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
! only for the nonlinear fields.
  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(4) = fifth
   if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

!Right Eigenvectors
  Rv(1,1) = one    
  Rv(2,1) = u - a*nx
  Rv(3,1) = v - a*ny
  Rv(4,1) = H - un*a

  Rv(1,2) = zero
  Rv(2,2) = mx
  Rv(3,2) = my
  Rv(4,2) = um

  Rv(1,3) = one
  Rv(2,3) = u
  Rv(3,3) = v 
  Rv(4,3) = half*(u*u+v*v)

  Rv(1,4) = one
  Rv(2,4) = u + a*nx
  Rv(3,4) = v + a*ny
  Rv(4,4) = H + un*a

!Dissipation Term
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*LdU(j)*Rv(i,j)
   end do
  end do

!Compute the flux.
  fL(1) = rhoL*unL
  fL(2) = rhoL*unL * uL + pL*nx
  fL(3) = rhoL*unL * vL + pL*ny
  fL(4) = rhoL*unL * HL

  fR(1) = rhoR*unR
  fR(2) = rhoR*unR * uR + pR*nx
  fR(3) = rhoR*unR * vR + pR*ny
  fR(4) = rhoR*unR * HR

  flux = half * (fL + fR - diss)
  wsn = half*(abs(un) + a)  !Normal max wave speed times half

 end subroutine roe
 
end module roe_scheme