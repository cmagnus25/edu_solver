# basic library
SET(BAS_src ${SRCBAS}/edu2d_euler_steady_basic_package_v0.f90
            ${SRCBAS}/gaussian_elimination.f90
            ${SRCBAS}/input_parameter.f90
            ${SRCBAS}/vector_operations.f90
)
ADD_LIBRARY(${BASICLIB} STATIC ${BAS_src})
