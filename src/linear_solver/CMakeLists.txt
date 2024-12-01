# linear solver library
SET(LINSOLV_src ${SRCLIN}/edu2d_euler_linear_solve_v0.f90
                ${SRCLIN}/frechet_derivative.f90
		${SRCLIN}/gcr_solver.f90
		${SRCLIN}/preconditioner.f90	
)
ADD_LIBRARY(${LINSOLVLIB} STATIC ${LINSOLV_src})

TARGET_LINK_LIBRARIES(${LINSOLVLIB} ${BASICLIB})
TARGET_LINK_LIBRARIES(${LINSOLVLIB} ${RESLIB})
