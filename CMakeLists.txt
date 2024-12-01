# CMake project file for FOO

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(EDU2D Fortran)

# Set the FOO version
SET(VERSION 1.2.3)

# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
#OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
#OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
# Locate and set parallelization libraries.  There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
#INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache. 
#INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
SET(MAINEXE edu2d)

# Define the library name
#SET(BARLIB bar)
SET(BASICLIB basic)
SET(GRADLIB gradient)
SET(LINSOLVLIB linear_solver)
SET(JACLIB jacobian)
#SET(RESLIB residual)
SET(SOLVLIB solver)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB ${CMAKE_SOURCE_DIR}/lib)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(SRCMAIN ${SRC})
SET(SRCBAS  ${SRC}/basic)
SET(SRCGRAD ${SRC}/gradient)
SET(SRCLIN  ${SRC}/linear_solver)
SET(SRCJAC  ${SRC}/jacobian)
#SET(SRCRES  ${SRC}/residual)
SET(SRCSOL  ${SRC}/solver)

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

# The source for the BAR library and have it placed in the lib folder
ADD_SUBDIRECTORY(${SRCBAS}  ${LIB}/basic)
ADD_SUBDIRECTORY(${SRCGRAD} ${LIB}/gradient)
ADD_SUBDIRECTORY(${SRCLIN}  ${LIB}/linear_solver)
ADD_SUBDIRECTORY(${SRCJAC}  ${LIB}/jacobian)
#DD_SUBDIRECTORY(${SRCRES}  ${LIB}/residual)
ADD_SUBDIRECTORY(${SRCSOL}  ${LIB}/solver)
# The source for the FOO binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRCMAIN} ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)
