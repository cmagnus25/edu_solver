# residual library
SET(RES_src ${SRCRES}/residual.f90
)
ADD_LIBRARY(${RESLIB} STATIC ${RES_src})
