# solver library
SET(SOLV_src ${SRCSOL}/edu2d_euler_steady_solver_v0.f90
)
ADD_LIBRARY(${SOLVLIB} STATIC ${SOLV_src})

TARGET_LINK_LIBRARIES(${SOLVLIB} ${BASICLIB})
TARGET_LINK_LIBRARIES(${SOLVLIB} ${GRADLIB})
TARGET_LINK_LIBRARIES(${SOLVLIB} ${JACLIB})
TARGET_LINK_LIBRARIES(${SOLVLIB} ${LINSOLVLIB})
TARGET_LINK_LIBRARIES(${SOLVLIB} ${RESLIB})
