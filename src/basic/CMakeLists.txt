# basic library
SET(BAS_src ${SRCBAS}/edu2d_euler_steady_basic_package_v0.f90
            ${SRCBAS}/gaussian_elimination.f90
)
ADD_LIBRARY(${BASICLIB} STATIC ${BAS_src})
