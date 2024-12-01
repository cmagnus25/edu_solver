# linear solver library
SET(LINSOLV_src ${SRCLIN}/edu2d_euler_linear_solve_v0.f90
)
ADD_LIBRARY(${LINSOLVLIB} STATIC ${LINSOLV_src})

TARGET_LINK_LIBRARIES(${LINSOLVLIB} ${BASICLIB})
