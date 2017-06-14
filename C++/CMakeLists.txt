# dmrosen 18-May-2017

# PROJECT CONFIGURATION
cmake_minimum_required(VERSION 3.1)
project(SESync)
set(CMAKE_CXX_STANDARD 11)  # We require C++ 11


### CMake Cache (build configuration) variables -- these are set interactively in the CMake GUI, and cached in CMakeCache ###

# Build type
SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
# Directory for built libraries
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib CACHE PATH "The directory in which to place libraries built by this project")
# Directory for built executables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin CACHE PATH "The directory in which to place executables built by this project")


# INCLUDE DIRECTORIES
set(SPECTRA_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/spectra/include CACHE PATH "Path to top-level directory of Spectra template library")
find_path(EIGEN_INCLUDE_DIR signature_of_eigen3_matrix_library HINTS ${CMAKE_SOURCE_DIR}/eigen DOC "Path to top-level directory of Eigen template library")
find_path(SUITESPARSE_INCLUDE_DIR SuiteSparse_config.h HINTS /usr/include/suitesparse)

# REQUIRED LIBRARIES
find_library(SUITESPARSECONFIG_LIB suitesparseconfig DOC "Suitesparseconfig library")
find_library(SPQR_LIB spqr DOC "SPQR library")
find_library(CHOLMOD_LIB cholmod DOC "Cholmod library")
set(SuiteSparse_LIBS ${SUITESPARSECONFIG_LIB} ${SPQR_LIB} ${CHOLMOD_LIB})



# Import the ROPTLIB library
add_subdirectory(ROPTLIB)

set(SESync_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SE-Sync/include)
set(SESync_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SE-Sync/src)

# Expose the include directories for this project
set(SESYNC_INCLUDE_DIRS ${SESync_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR} ${SPECTRA_INCLUDE_DIR} ${SUITESPARSE_INCLUDE_DIR} ${ROPTLIB_INCLUDE_DIRECTORIES} CACHE INTERNAL "")
include_directories(${SESync_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR} ${SPECTRA_INCLUDE_DIR} ${SUITESPARSE_INCLUDE_DIR} ${ROPTLIB_INCLUDE_DIRECTORIES})


# Get the set of SE-Sync header and source files
set(SESync_HDRS
${SESync_INCLUDE_DIR}/RelativePoseMeasurement.h
${SESync_INCLUDE_DIR}/SESync_utils.h
${SESync_INCLUDE_DIR}/SESyncProblem.h
${SESync_INCLUDE_DIR}/SESync.h
)

set(SESync_SRCS
${SESync_SOURCE_DIR}/SESync_utils.cpp
${SESync_SOURCE_DIR}/SESyncProblem.cpp
${SESync_SOURCE_DIR}/SESync.cpp
${SESync_SOURCE_DIR}/main.cpp
)

# Build the SE-Sync library
add_library(SESync ${SESync_HDRS} ${SESync_SRCS})
target_link_libraries(SESync roptlib ${SuiteSparse_LIBS})

# Build the example executable
add_executable(se-sync ${SESync_SOURCE_DIR}/main.cpp)
target_link_libraries(se-sync SESync ${SuiteSparse_LIBS})

