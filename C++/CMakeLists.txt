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

# BUILD CONFIGURATIONS
option(CMAKE_VERBOSE_MAKEFILE "Generate verbose makefiles?" OFF)

# Add the .cmake files that ship with Eigen3 to the CMake module path (useful for finding other stuff)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/eigen/cmake" CACHE STRING "The CMake module path used for this project")


# SESYNC DEPENDENCIES

# FIND EIGEN3
find_package(Eigen3 3.3 REQUIRED)
if(EIGEN3_FOUND)
message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION_STRING})")
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
message(STATUS "Eigen library not found!")
endif()

# FIND ADDITIONAL LIBRARIES
# These next operations make use of the .cmake files shipped with Eigen3 
find_package(SPQR REQUIRED)
find_package(Cholmod REQUIRED)
find_package(BLAS REQUIRED)

# FIND SPECTRA
# Spectra doesn't ship with a .cmake file, but it's header-only, so it suffices to simply determine it's top-level directory
set(SPECTRA_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/spectra/include CACHE PATH "Path to top-level directory of Spectra template library")


# IMPORT ROPTLIB library
add_subdirectory(ROPTLIB)

set(SESync_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SE-Sync/include)
set(SESync_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SE-Sync/src)
set(SESync_EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/examples)

# Expose the include directories for this project
set(SESYNC_INCLUDES ${SESync_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${SPQR_INCLUDES} ${CHOLMOD_INCLUDES} ${SPECTRA_INCLUDE_DIR}  ${ROPTLIB_INCLUDE_DIRECTORIES} CACHE INTERNAL "")



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
)

# Build the SE-Sync library
add_library(SESync SHARED ${SESync_HDRS} ${SESync_SRCS})
target_include_directories(SESync PUBLIC ${SESYNC_INCLUDES})
target_link_libraries(SESync roptlib ${BLAS_LIBRARIES} ${CHOLMOD_LIBRARIES} ${SPQR_LIBRARIES} ${M} ${LAPACK})

# Build the example executable
add_executable(se-sync ${SESync_EXAMPLES_DIR}/main.cpp)
target_include_directories(se-sync PUBLIC ${SESYNC_INCLUDES})
target_link_libraries(se-sync SESync ${BLAS_LIBRARIES} ${CHOLMOD_LIBRARIES} ${SPQR_LIBRARIES} ${M} ${LAPACK})

# Export the SESync and ROPTLIB libraries to make them available to other CMake projects
export(TARGETS roptlib SESync FILE SESync.cmake)
