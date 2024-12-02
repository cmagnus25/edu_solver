!********************************************************************************
! Educationally-Designed Unstructured 2D (EDU2D) Code
!
!  ---------------- EDU2D-Euler-IMPLICIT
!
!  - Input parameter module
!
! This module defines input parameters, and set the default values.
! These parameters are specified in the file named 'input.nml', and read
! in the main program by the subroutine "read_nml_input_parameters".
!
!
!        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!
! the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!
! 
! This F90 program is written and made available for an educational purpose.
!
! This file may be updated in future.
!
! Katate Masatsuka http://www.cfdbooks.com
!********************************************************************************

 module input_parameter

  use edu2d_constants , only : p2, one, zero

  implicit none

 !This module contains the following subroutine:
 ! -read_nml_input_parameters ! read the input file

 !Make all input parameters andn subroutines available to other modules.
  public

!--------------------------------------------------------------------------------
! Set file names

! Input files
  character(80) :: datafile_grid_in  = 'bump.grid'
  character(80) :: datafile_bcmap_in = 'bump.bcmap'

! Output file
  character(80) :: datafile_tec      = 'bump_solution_tecplot.dat'

!--------------------------------------------------------------------------------
! Input Parameters

  real(p2) :: M_inf  = 0.3_p2     ! Free stream Mach number
  real(p2) :: gamma  = 1.4_p2     ! Ratio of specific heats

  character(80) :: iteration_method = "implicit_gcr" ! Solution method: "explicit" or "implicit"
  character(80) :: smooth_method    = "sgs"      ! gs or sgs, smoothing method
	 
!  For explicit scheme (2-staege Runge-Kutta)

  real(p2) :: CFLexp = 0.99_p2    ! CFL for expicit method (RK2); must be small

!  For implicit scheme (exact first-order Jacobian with a pseudo time term)
!  Note: Explore various combinations of these parameters.
!        E.g., Increase 'sweeps' to solve the linear system better.

  real(p2) :: CFL1 = 1.0e+1_p2  ! Initial CFL for implicit method
  real(p2) :: CFL2 = 1.0e+5_p2  !   Final CFL for implicit method
  real(p2) :: CFL_ramp_steps = 10         ! Number of iterations to reach from CFL1 to CFL2
  integer  :: sweeps = 30        ! Number of linear GS sweeps for implicit method

  real(p2) :: tolerance = 1.0e-15_p2 ! Residual tolerance for steady computations
  real(p2) :: tolerance_linear = 0.5e+0_p2  ! Residual tolerance for linear system
  real(p2) :: tolerance_gcr = 1.0e-1_p2  !
  integer  :: max_iterations = 20         ! Max number of iterations
  integer  :: max_projection_gcr = 30         ! Max projections for GCR (larger->more expensive/memory)
	  
!  Sorry, but only the Roe flux is implemented in this code.

  character(80) :: inviscid_flux = "roe"      ! "roe" - Roe flux only. You can add others.
  character(80) :: inviscid_jac  = "roe"      ! "roe" - Roe flux Jacobian only. You can add others.

  character(80) :: gradient_type     = "linear" ! or "quadratic2" for a quadratic LSQ.
  character(80) :: gradient_weight   = "none"   ! or "inverse_distance"
  real(p2)      :: gradient_weight_p = zero    ! or any other real value

! End of Default input values
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
! Group them into "input_parameters":

 namelist / input_parameters / &
  datafile_grid_in    , &
  datafile_bcmap_in   , &
  datafile_tec        , &
  M_inf               , &
  gamma               , &
  iteration_method    , &
  smooth_method       , &
  CFLexp              , &
  CFL1                , & 
  CFL2                , &
  CFL_ramp_steps      , &
  sweeps              , &
  tolerance           , &
  tolerance_linear    , &
  tolerance_gcr       , &
  max_iterations      , &
  max_projection_gcr  , &
  inviscid_flux       , &
  inviscid_jac        , &
  gradient_type       , &
  gradient_weight     , &
  gradient_weight_p

  !Note: These variables defined above are available within the entire
  !      module, and so within all subroutines contained below.

 contains

!*****************************************************************************
!
! This subroutine reads input_parameters specified in the input file:
!
!    file name = namelist_file
!
! prints the content on screen.
!
! In the main program, we set: namelist_file = "input.nml"
!
!*****************************************************************************

  subroutine read_nml_input_parameters(namelist_file)

  implicit none
  character(9), intent(in) :: namelist_file
  integer :: os

  write(*,*)
  write(*,*) "-------------------------------------------------------"
  write(*,*) " Reading the input file: input.nml..... "
  write(*,*)

  open(unit=10,file=trim(namelist_file),form='formatted',status='old',iostat=os)
  read(unit=10,nml=input_parameters)

  write(*,*)
  write(*,*) " List of given namelist variables and their values"
  write(*,*)

  write(*,nml=input_parameters) ! Print the namelist variables on screen.
  close(10)

  write(*,*)
  write(*,*) " End of Reading the input file: input.nml..... "
  write(*,*) "-------------------------------------------------------"
  write(*,*)

  end subroutine read_nml_input_parameters
!*****************************************************************************

 end module input_parameter
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!  End of input parameter module
