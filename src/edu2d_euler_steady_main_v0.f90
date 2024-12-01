!********************************************************************************
!* Educationally-Designed Unstructured 2D (EDU2D) Code
!*
!*
!*                     --- This is EDU2D-Euler-Implicit ---
!*
!*
!* EDU2D-Euler-Implicit: An Euler code with
!*
!*    - Node-centered finite-volume discretization
!*    - Fully implicit solver with the exact first-order Jacobian
!*
!*
!*             written especially for steady state computations
!*
!* - Node-centered second-order finite-volume method for unstructured grids
!* - The Roe flux with an entropy fix
!* - Exact Roe flux Jacobian for implicit formulation
!* - Gradient reconstruction by unweighted least-squares method
!* - Explicit: 2-Stage Runge-Kutta local time-stepping towards steady state
!* - Implicit: Defect correction (DC) method: linearization based on 1st-order scheme
!* - Sequential Gauss-Seidel relaxation for the linear system
!* - All quantities are nondimensionalized; velocity and pressure are
!*   nondimensionalized based on the free stream speed of sound
!*   (see Section 4.8.3 in I do like CFD, VOL.1).
!* - Boudnary flux quadrature is implemented specially for triangular grids.
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!*
!* This is Version 0 (July 2015)
!*
!* ------------------------------------------------------------------------------
!* Files: There are 5 files.
!*
!* ------------------------------------------
!* - Main driver program file   : This reads grid and BC files, and call dummy NC/CC programs.
!*
!*     edu2d_euler_implicit_main.f90 (This file), which contains a main driver program
!*      -- edu2d_euler_implicit : Main driver code, which calls an Euler solver
!*
!* ------------------------------------------
!* - Basic EDU2D package file   : Arranged for a 2D implicit Euler code
!*
!*     edu2d_basic_package_euler_implct.f90, which contains the following modules.
!*      -- edu2d_constants      : Numerical values defined
!*      -- edu2d_grid_data_type : Grid data types defined
!*      -- edu2d_main_data      : Main grid data and parameters declared
!*      -- edu2d_grid_data      : Read/construct/check grid data
!*
!* ------------------------------------------
!* - Euler solver file   : This computes a solution to the shock diffraction problem.
!*
!*     edu2d_euler_implct_solver.f90, which contains a 2D Euler solver with RK2.
!*      -- edu2d_euler_implct_solver : Node-centered Explicit Euler solver with RK2
!*
!* ------------------------------------------
!* - Jacobian file   : This computes a solution to the shock diffraction problem.
!*
!*     edu2d_euler_jacobian.f90, which contains a 2D Euler solver with RK2.
!*      -- edu2d_euler_jacobian : Node-centered Explicit Euler solver with RK2
!*
!* ------------------------------------------
!* - Linear solver file   : This computes a solution to the shock diffraction problem.
!*
!*     edu2d_euler_linear_solve.f90, which contains a 2D Euler solver with RK2.
!*      -- edu2d_euler_linear_solve : Node-centered Explicit Euler solver with RK2
!*
!* ------------------------------------------------------------------------------
!* Notes:
!*
!*  The purpose of this code is to give a beginner an opportunity to learn how to
!*  write an unstructured CFD code. Hence, the focus here is on the simplicity.
!*  The code is not optimized for efficiency.
!*
!*  If the code is not simple enough to understand, please send questions to Hiro
!*  at sunmasen(at)hotmail.com. I'll greatly appreciate it and revise the code.
!*
!*  If the code helps you understand how to write your own code that is more
!*  efficient and has more features, it'll have served its purpose.
!*
!* ------------------------------------------------------------------------------
!* Examples of additional features you might want to add.
!*
!*  1. More boundary conditions (periodic, symmetry, suction, etc.)
!*  2. Other reconstruction     (Van Leer's kappa scheme, etc.)
!*  3. Limiters                 (Venkat/Barth limiter,etc.)
!*  4. Other flux functions     (HLL, LDFSS, AUSM, etc.)
!*  5. Local-preconditioning    (low-Mach number accuracy and stiffness removal)
!*  6. More output              (convergence history, etc.)
!*  7. Parallelization          (large-scale problems)
!*  8. Grid adaptation          (h-refinement, steady or unsteady)
!*  9. Adjoint capability       (aerodynamic design, adaptation, etc.)
!* 10. Moving grid              (sliding mesh, overset grid, etc.)
!* 11. Multigrid                (grid-independent convergence)
!* ------------------------------------------------------------------------------
!*
!* Katate Masatsuka, July 2015. http://www.cfdbooks.com
!********************************************************************************

!********************************************************************************
!* Main program: Node-centered finite-volume Euler code
!*
!* This code computes steady solutions of the Euler equations.
!*
!* Input -------------------------------------------------------
!*
!*   bump.grid  = grid file containing boundary information
!*   bump.bcmap = file that specifies boundary conditions
!*
!*   (Note: See the subroutine "read_grid", which is in this file,
!*          for the format of these files.)
!*
!* Output ------------------------------------------------------
!*
!*  "bump_solution_tecplot.dat" = Tecplot file containing the grid and the solution.
!*
!*
!*  NOTE: The variables are nondimensionalized values (compressible scale),
!*           rho=rho/rho_inf, u=u/a_inf, v=v/a_inf, rho=p/(rho_inf*a_inf^2)
!*
!*  NOTE: Many variables are passed to subroutines via USE statement.
!*        Each module contains specific data, and they are accessed by USE.
!*
!*
!* Katate Masatsuka, July 2015. http://www.cfdbooks.com
!********************************************************************************
 program ossan_euler2d

 use edu2d_constants    , only : p2, zero
 use edu2d_grid_data    , only : read_grid, construct_grid_data, check_grid_data
 use edu2d_my_main_data , only : nq, M_inf, gamma, iteration_method,         &
                                 CFLexp, CFL1, CFL2, CFL_ramp_steps, sweeps, &
                                 tolerance, tolerance_linear, max_iterations,&    
                                 inviscid_flux, inviscid_jac,                &
                                 jac, nnodes, node, gradient_type, gradient_weight, gradient_weight_p

 use edu2d_euler_implct_solver, only : euler_solver_main
 use gradients_lsq, only : compute_lsq_coeff_nc, check_lsq_coeff_nc
 
 implicit none

!Inout data files
 character(80) :: datafile_grid_in  !Grid file
 character(80) :: datafile_bcmap_in !Boundary condition specification file

!Output data file
 character(80) :: datafile_tec      !Tecplot file for viewing the result

!Local variables
 integer       :: i

  write(*,*) "***************************************************************"
  write(*,*) " Starting EDU2D-Euler-Steady "
  write(*,*) "***************************************************************"

!--------------------------------------------------------------------------------
! Set file names

! Input files
  datafile_grid_in  = 'bump.grid'
  datafile_bcmap_in = 'bump.bcmap'

! Output file
  datafile_tec      = 'bump_solution_tecplot.dat'

!--------------------------------------------------------------------------------
! Input Parameters

            M_inf  = 0.3_p2     ! Free stream Mach number
             gamma = 1.4_p2     ! Ratio of specific heats

  iteration_method = "implicit" ! Solution method: "explicit" or "implicit"

!  For explicit scheme (2-staege Runge-Kutta)

            CFLexp = 0.99_p2    ! CFL for expicit method (RK2); must be small

!  For implicit scheme (exact first-order Jacobian with a pseudo time term)
!  Note: Explore various combinations of these parameters.
!        E.g., Increase 'sweeps' to solve the linear system better.

              CFL1 = 1.0e+00_p2 ! Initial CFL for implicit method
              CFL2 = 1.0e+5_p2  !   Final CFL for implicit method
    CFL_ramp_steps = 120        ! Number of iterations to reach from CFL1 to CFL2
            sweeps = 15         ! Number of linear GS sweeps for implicit method

         tolerance = 1.0e-15_p2 ! Residual tolerance for steady computations
  tolerance_linear = 1.0e-1_p2 ! Residual tolerance for linear system
    max_iterations = 1000      ! Max number of iterations

!  Sorry, but only the Roe flux is implemented in this code.

     inviscid_flux = "roe"      ! "roe" - Roe flux only. You can add others.
     inviscid_jac  = "roe"      ! "roe" - Roe flux Jacobian only. You can add others.

                  nq = 4        ! The number of equtaions/variables in the target equtaion.
                                ! This is 4 for 2D Euler equations.

   gradient_type     = "linear" ! or "quadratic2 for a quadratic LSQ.
   gradient_weight   = "none"   ! or "inverse_distance"
   gradient_weight_p =  zero    ! or any other real value

!-----------------------------------------------------------------------------
! Note: The implicit solver becomes Newton's method for first-order scheme.
!       You can try it with the following setting:
!
!             CFL1 = 1.0e+15_p2
!             CFL2 = 1.0e+15_p2
!           sweeps = 250
!    gradient_type = "none"
!
!-----------------------------------------------------------------------------

! Print out the input parameters just for a record.

  write(*,*)
  write(*,*) " Summary of Input Parameteres:"
  write(*,*)
  write(*,*) "                  M_inf = ", M_inf
  write(*,*) "                  gamma = ", gamma
  write(*,*)
  write(*,*) "       iteration_method = ", trim(iteration_method)
  write(*,*)
  write(*,*) "                 CFLexp = ", CFLexp !CFL for explicit method
  write(*,*)
  write(*,*) "                   CFL1 = ", CFL1
  write(*,*) "                   CFL2 = ", CFL2
  write(*,*) "         CFL_ramp_steps = ", CFL_ramp_steps
  write(*,*) "                 sweeps = ", sweeps
  write(*,*)
  write(*,*) "              tolerance = ", tolerance
  write(*,*) "       tolerance_linear = ", tolerance_linear
  write(*,*) "         max_iterations = ", max_iterations
  write(*,*)
  write(*,*) "          inviscid_flux = ", trim(inviscid_flux)
  write(*,*) "          inviscid_jac  = ", trim(inviscid_jac)
  write(*,*)
  write(*,*) "          gradient_type = ", trim(gradient_type)
  write(*,*) "        gradient_weight = ", trim(gradient_weight)
  write(*,*) "      gradient_weight_p = ", gradient_weight_p
  write(*,*)

!--------------------------------------------------------------------------------
! Read a grid, solve the Euler equations, and write the output datafile.
!

! (1) Read grid files
      call read_grid(datafile_grid_in, datafile_bcmap_in)

! (2) Construct grid data
      call construct_grid_data

     !Allocate additional arrays needed for the Euler solver.
      do i = 1, nnodes
       allocate( node(i)%u(    nq  ) )
       allocate( node(i)%du(   nq  ) )
       allocate( node(i)%w(    nq  ) )
       allocate( node(i)%gradw(nq,2) ) !<- 2: x and y components.
       allocate( node(i)%res(  nq  ) )
      end do

     !For implicit method, allocate the arrays that store a Jacobian matrix
      if (trim(iteration_method) == "implicit") then
         allocate(jac(nnodes))
       do i = 1, nnodes
         allocate( jac(i)%diag(nq,nq)                 ) ! Diagonal block
		 allocate( jac(i)%diag_inverse(nq,nq)         ) ! Inverse Diagonal block
         allocate( jac(i)%off( nq,nq, node(i)%nnghbrs)) ! Off-diagonal block
       end do
      endif

! (3) Check the grid data (It is always good to check them before use!)
      call check_grid_data

! (4) Compute, store, and check LSQ coefficients
      call compute_lsq_coeff_nc
      call check_lsq_coeff_nc

! (5) Set initial solution
      call set_initial_solution

! (6) Compute the numerical solution
      call euler_solver_main

! (7) Write out the tecplot data file (Solutions at nodes)
      call write_tecplot_file(datafile_tec)

!--------------------------------------------------------------------------------

  write(*,*) "Successfully completed. Stop."

 stop

 contains

!********************************************************************************
!* Initial solution:
!*
!* Set initial solution here (typically free stream condition).
!********************************************************************************
 subroutine set_initial_solution

 use edu2d_constants   , only : zero, one
 use edu2d_my_main_data, only : nnodes, node, gamma, M_inf, rho_inf, u_inf, v_inf, p_inf
 use edu2d_euler_implct_solver, only : w2u

 implicit none

!Local variables
 integer  :: i

 write(*,*)
 write(*,*) "Setting up initial solution (uniform stream)..."
 write(*,*)

! Uniform stream values (nondimensionalized variables)

   rho_inf = one
     u_inf = M_inf
     v_inf = zero
     p_inf = one/gamma

! Specify the uniform stream values at all nodes.

  nodes : do i = 1, nnodes

   node(i)%w = (/ rho_inf, u_inf, v_inf, p_inf /) ! Primitive variables
   node(i)%u = w2u(node(i)%w)                     ! Conservative variables

  end do nodes

 end subroutine set_initial_solution

!********************************************************************************
!* Write a tecplot file: grid and solution
!********************************************************************************
 subroutine write_tecplot_file(datafile_tec)

 use edu2d_my_main_data, only : nnodes, node, elm, nelms, gamma

 implicit none
 character(80), intent(in) :: datafile_tec
 integer :: i, k, os
!--------------------------------------------------------------------------------
 open(unit=1, file=datafile_tec, status="unknown", iostat=os)

 write(1,*) 'title = "grid"'

 write(1,'(a80)') 'variables = "x","y","rho","u","v","p","Mach"'

 write(1,*) 'zone n=',nnodes,' e =', nelms,' et=quadrilateral, f=fepoint'

!--------------------------------------------------------------------------------
! Nodal quantities: x, y, rho, u, v, p, Mach number

   do i = 1, nnodes
    write(1,*) node(i)%x, node(i)%y, (node(i)%w(k),k=1,4), &
               sqrt( (node(i)%w(2)**2+node(i)%w(3)**2)/(gamma*node(i)%w(4)/node(i)%w(1)) )
   end do

!--------------------------------------------------------------------------------
! Both quad and tria elements in quad format:

 do i = 1, nelms

  if (elm(i)%nvtx == 3) then

   write(1,*) elm(i)%vtx(1), elm(i)%vtx(2), elm(i)%vtx(3), elm(i)%vtx(3)

  elseif (elm(i)%nvtx == 4) then

   write(1,*) elm(i)%vtx(1), elm(i)%vtx(2), elm(i)%vtx(3), elm(i)%vtx(4)

  else

   !Impossible
   write(*,*) " Error in elm%vtx data... Stop..: elm(i)%nvtx=",elm(i)%nvtx
   stop

  endif

 end do

!--------------------------------------------------------------------------------
 close(1)
 end subroutine write_tecplot_file
!********************************************************************************

 end program ossan_euler2d




