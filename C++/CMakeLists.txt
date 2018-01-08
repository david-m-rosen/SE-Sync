# dmrosen 18-May-2017

cmake_minimum_required(VERSION 3.1)

# PROJECT CONFIGURATION
project(SESync C CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # We require C++ 14

# Enable faster instruction sets (SIMD/AVX)
add_compile_options(-march=native)

# Set build type to 'RelWithDebInfo' if one was not specified by the user
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
  message(STATUS "Setting build type to ${CMAKE_BUILD_TYPE}, as none was specified\n")
else()
  message(STATUS "Building in ${CMAKE_BUILD_TYPE} mode\n")
endif()



# Directory for built libraries
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib CACHE PATH "The directory in which to place the SE-Sync library built by this project")
# Directory for built executables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin CACHE PATH "The directory in which to place executables built by this project")

# BUILD CONFIGURATIONS
option(CMAKE_VERBOSE_MAKEFILE "Generate verbose makefiles?" OFF)

# Add the .cmake files that ship with Eigen3 to the CMake module path (useful for finding other stuff)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/eigen/cmake" CACHE STRING "The CMake module path used for this project")


# SESYNC DEPENDENCIES

# FIND EIGEN3
set( ENV{EIGEN3_ROOT_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/eigen)
find_package(Eigen3 3.3.3 REQUIRED)
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


#Find Optimization library
set(OPTIMIZATION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Optimization/ CACHE PATH "Path to top-level directory of Optimization (the one containing CMakeLists.txt)")
# Run this project's CMakeLists file
add_subdirectory(${OPTIMIZATION_DIR})
#Import dependencies from this library
find_library(OptLib NAMES Optimization)

# PERFORMANCE IMPROVEMENTS
# Uncomment the following 'find_package' line to turn on support for multithreading (if available)
# Check if OpenMP support is available (for parallelization of Eigen operations)
find_package(OpenMP)
if(OPENMP_FOUND)
message(STATUS "\nFound OpenMP! Turning on support for parallelization")
endif()


# FIND SPECTRA
# Spectra doesn't ship with a .cmake file, but it's header-only, so it suffices to simply determine it's top-level directory
set(SPECTRA_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/spectra/include CACHE PATH "Path to top-level directory of Spectra template library")

set(SESync_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SE-Sync/include)
set(SESync_HDR_DIR ${SESync_INCLUDE_DIR}/SESync)
set(SESync_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SE-Sync/src)
set(SESync_EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/examples)

# Expose the include directories for this project
set(SESYNC_INCLUDES ${SESync_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${SPQR_INCLUDES} ${CHOLMOD_INCLUDES} ${SPECTRA_INCLUDE_DIR} CACHE INTERNAL "")


# SE-SYNC PROJECT

# Get the set of SE-Sync header and source files
set(SESync_HDRS
${SESync_HDR_DIR}/StiefelProduct.h
${SESync_HDR_DIR}/RelativePoseMeasurement.h
${SESync_HDR_DIR}/SESync_types.h
${SESync_HDR_DIR}/SESync_utils.h
${SESync_HDR_DIR}/SESyncProblem.h
${SESync_HDR_DIR}/SESync.h
)

set(SESync_SRCS
${SESync_SOURCE_DIR}/StiefelProduct.cpp
${SESync_SOURCE_DIR}/SESync_utils.cpp
${SESync_SOURCE_DIR}/SESyncProblem.cpp
${SESync_SOURCE_DIR}/SESync.cpp
)

# Build the SE-Sync library
add_library(SESync SHARED ${SESync_HDRS} ${SESync_SRCS})
target_include_directories(SESync PUBLIC ${SESYNC_INCLUDES})
target_link_libraries(SESync Optimization ${BLAS_LIBRARIES} ${CHOLMOD_LIBRARIES} ${SPQR_LIBRARIES} ${M} ${LAPACK})

# Uncomment these lines to enable support for OpenMP
if(OPENMP_FOUND)
# Add additional compilation flags to enable OpenMP support
set_target_properties(SESync PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
target_link_libraries(SESync gomp)
endif()

# Export the SESync library to make it available to other CMake projects
export(TARGETS SESync FILE SESync.cmake)

# BUILD EXAMPLE DRIVER
add_subdirectory(examples)

