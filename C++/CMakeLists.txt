# dmrosen 18-May-2017

cmake_minimum_required(VERSION 3.1)

# PROJECT CONFIGURATION
project(SESync LANGUAGES C CXX VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # We require C++ 14


# Set build type to 'RelWithDebInfo' if one was not specified by the user
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS RelWithDebInfo Release Debug MinSizeRel)
  message(STATUS "Setting build type to ${CMAKE_BUILD_TYPE}, as none was specified\n")
else()
  message(STATUS "Building in ${CMAKE_BUILD_TYPE} mode\n")
endif()

# Directory for built libraries
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib CACHE PATH "The directory in which to place the SE-Sync library built by this project")
# Directory for built executables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin CACHE PATH "The directory in which to place executables built by this project")

# BUILD CONFIGURATIONS
# Enable faster instruction sets (SIMD/AVX)
set(ENABLE_VECTORIZATION OFF CACHE BOOL "Enable vectorized instruction sets (SIMD/AVX)? [disabled by default]")
# Enable OpenMP (if available)
set(ENABLE_OPENMP ON CACHE BOOL "Enable OpenMP (if available)")
# Enable code profiling using gperftools
set(ENABLE_PROFILING OFF CACHE BOOL "Enable code profiling using gperftools")
# Enable visualization module.
set(ENABLE_VISUALIZATION OFF CACHE BOOL "Enable visualization module.")
# Build Python bindings
set(BUILD_PYTHON_BINDINGS OFF CACHE BOOL "Build Python bindings.")

# Add the .cmake files that ship with Eigen3 to the CMake module path (useful for finding other stuff)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Eigen/cmake" CACHE STRING "The CMake module path used for this project")

if(${ENABLE_VECTORIZATION})
message(STATUS "Enabling SIMD/AVX instruction sets")
add_definitions(-march=native)
endif()

if(${ENABLE_PROFILING})
message(STATUS "Enabling code profiling using Google Performance Tools")
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -lprofiler)
endif()

if(${ENABLE_VISUALIZATION})
  find_package(Pangolin REQUIRED)
endif()

if(${BUILD_PYTHON_BINDINGS})
message(STATUS "Building Python bindings")
endif()


message(STATUS "")


# SESYNC DEPENDENCIES

# FIND EIGEN3
set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Eigen)
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")


# FIND ADDITIONAL LIBRARIES
# These next operations make use of the .cmake files shipped with Eigen3
find_package(SPQR REQUIRED)
find_package(Cholmod REQUIRED)
find_package(BLAS REQUIRED)


# Find Optimization library
set(OPTIMIZATION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Optimization/ CACHE PATH "Path to top-level directory of Optimization (the one containing CMakeLists.txt)")

# Run the Optimization package's CMake file
add_subdirectory(${OPTIMIZATION_DIR})

# PERFORMANCE IMPROVEMENTS
if(${ENABLE_OPENMP})
find_package(OpenMP)
if(OPENMP_FOUND)
message(STATUS "Found OpenMP! Turning on support for parallelization\n")
endif()
endif()


# FIND SPECTRA
# Spectra doesn't ship with a .cmake file, but it's header-only, so it suffices to simply determine it's top-level directory
set(SPECTRA_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/include CACHE PATH "Path to top-level directory of Spectra template library")
message(STATUS "Spectra include directory: ${SPECTRA_INCLUDE_DIR}\n")

set(SESync_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SE-Sync/include)
set(SESync_HDR_DIR ${SESync_INCLUDE_DIR}/SESync)
set(SESync_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SE-Sync/src)
set(SESync_EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/examples)

# Set the include directories for this project
# Spectra is only used internally by the SE-Sync library
set(SESync_PRIVATE_INCLUDES ${SPECTRA_INCLUDE_DIR} CACHE INTERNAL "")
# The SE-Sync headers and and Eigen 3, SPQR, and Cholmod are all referenced by the header files of the SE-Sync library, hence must be PUBLICLY included (i.e. clients using the SE-Sync headers must also include these headers)
set(SESync_INCLUDES ${SESync_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${SPQR_INCLUDES} ${CHOLMOD_INCLUDES} CACHE INTERNAL "")


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
add_library(${PROJECT_NAME} SHARED ${SESync_HDRS} ${SESync_SRCS})
target_include_directories(${PROJECT_NAME} PRIVATE ${SESync_PRIVATE_INCLUDES})
target_include_directories(${PROJECT_NAME} PUBLIC ${SESync_INCLUDES})
target_link_libraries(${PROJECT_NAME} Optimization ${BLAS_LIBRARIES} ${CHOLMOD_LIBRARIES} ${SPQR_LIBRARIES} ${M} ${LAPACK})

if(OPENMP_FOUND)
# Add additional compilation flags to enable OpenMP support
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "-fopenmp")
endif()


# SE-SYNC VISUALIZATION PROJECT

# If visualization is enabled, build library.
if(${ENABLE_VISUALIZATION})
  set(SESyncViz_HDRS ${SESync_HDR_DIR}/SESyncVisualizer.h)
  set(SESyncViz_SRCS ${SESync_SOURCE_DIR}/SESyncVisualizer.cpp)

  add_library(SESyncViz SHARED ${SESyncViz_HDRS} ${SESyncViz_SRCS})
  target_include_directories(SESyncViz PUBLIC ${SESync_INCLUDE_DIR})
  target_link_libraries(SESyncViz ${PROJECT_NAME} ${Pangolin_LIBRARIES})
endif()

# PYTHON BINDINGS
if(${BUILD_PYTHON_BINDINGS})
  
  # Find the appropriate Python3 interpreter
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  
  # Find pybind11
  find_package(pybind11 REQUIRED)
  
  # Determine the location of this user's site-packages directory
  execute_process ( COMMAND ${Python3_EXECUTABLE} -m site --user-site OUTPUT_VARIABLE Python3_USER_SITE_LIB OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Add the Python module to this project
  pybind11_add_module(PySESync MODULE ${SESync_SOURCE_DIR}/PySESync.cpp)
  target_link_libraries(PySESync PRIVATE ${PROJECT_NAME})
  
  # Build this Python library in the user's local site-packages directory (as reported by Python itself), so that it can be automatically be found and imported by the Python interpreter
  message(STATUS "Building SE-Sync Python library in ${Python3_USER_SITE_LIB}")
  set_target_properties(PySESync PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${Python3_USER_SITE_LIB})
endif()



# BUILD EXAMPLE DRIVER
add_subdirectory(examples)


# EXPORT SE-SYNC LIBRARY

# Add add entry for this project into CMake's package registry, so that this project can be found by other CMake projects
export(PACKAGE ${PROJECT_NAME})
# Create a configuration file for this project, so that it can be imported by other CMake projects
export(TARGETS ${PROJECT_NAME} Optimization FILE ${PROJECT_NAME}Config.cmake)


# EXPORT SE-SYNC VISUALIZATION LIBRARY

if(${ENABLE_VISUALIZATION})
# Add add entry for this project into CMake's package registry, so that this
# project can be found by other CMake projects
export(PACKAGE SESyncViz)
# Create a configuration file for this project, so that it can be imported by
# other CMake projects
export(TARGETS SESyncViz FILE SESyncVizConfig.cmake)
endif()
