# jacobian library
SET(JAC_src ${SRCJAC}/edu2d_euler_jacobian_v1.f90
)
ADD_LIBRARY(${JACLIB} STATIC ${JAC_src})

TARGET_LINK_LIBRARIES(${JACLIB} ${BASICLIB})
