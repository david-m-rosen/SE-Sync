# dmrosen 17-May-2017:  This CMakeLists file is based upon the original ROPTLIB Makefile

# PROJECT CONFIGURATION
cmake_minimum_required(VERSION 3.1)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")

project(ROPTLIB C CXX)
set(CMAKE_CXX_STANDARD 11)  # We require C++ 11
add_compile_options(-march=native) # Enable faster instruction sets

# Build type

# Directory for built libraries
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib CACHE PATH "The directory in which to place libraries built by this project")
# Directory for built executables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin CACHE PATH "The directory in which to place executables built by this project")

# BUILD CONFIGURATIONS
option(CMAKE_VERBOSE_MAKEFILE "Generate verbose makefiles?" OFF)

set(CODE_PROFILING OFF CACHE BOOL "Turn on code profiling?")
if(${CODE_PROFILING})
message(STATUS "Turning on code profiling for ROPTLIB")
endif()


# ROPTLIB requires lapack, blas, and m libraries
find_library(LAPACK lapack DOC "LAPACK library to link against")
find_library(BLAS blas DOC "BLAS library to link against")
find_library(M m "M library to link against")

# Directories of ROPTLIB header files
set(ROPTLIB_INCLUDE_DIRECTORIES 
${CMAKE_CURRENT_SOURCE_DIR}
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds 
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/CpxNStQOrth 
#${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/ElasticShape 
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/EucPositive 
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/Euclidean 
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/Grassmann
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/L2Sphere
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/LowRank
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/Oblique
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/OrthGroup
#${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/PreShapeCurves
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/SPDManifold
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/Sphere
${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/Stiefel
${CMAKE_CURRENT_SOURCE_DIR}/Others
${CMAKE_CURRENT_SOURCE_DIR}/Others/SparseBLAS
${CMAKE_CURRENT_SOURCE_DIR}/Problems
${CMAKE_CURRENT_SOURCE_DIR}/Solvers
${CMAKE_CURRENT_SOURCE_DIR}/cwrapper
${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/blas
${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/lapack
CACHE INTERNAL ""
)

# COMPILE THE ROPTLIB LIBRARY
add_library(roptlib SHARED
Manifolds/Element.cpp Manifolds/LinearOPE.cpp Manifolds/Manifold.cpp Manifolds/ProductElement.cpp Manifolds/ProductManifold.cpp Manifolds/SharedSpace.cpp Manifolds/SmartSpace.cpp 
Manifolds/CpxNStQOrth/CSOVariable.cpp Manifolds/CpxNStQOrth/CSOVector.cpp Manifolds/CpxNStQOrth/CpxNStQOrth.cpp 
Manifolds/ElasticShape/ElasticShape.cpp Manifolds/ElasticShape/ShapeVariable.cpp Manifolds/ElasticShape/ShapeVector.cpp 
Manifolds/EucPositive/EucPosVariable.cpp Manifolds/EucPositive/EucPosVector.cpp Manifolds/EucPositive/EucPositive.cpp 
Manifolds/Euclidean/EucVariable.cpp Manifolds/Euclidean/EucVector.cpp Manifolds/Euclidean/Euclidean.cpp 
Manifolds/Grassmann/GrassVariable.cpp Manifolds/Grassmann/GrassVector.cpp Manifolds/Grassmann/Grassmann.cpp 
Manifolds/L2Sphere/L2Sphere.cpp Manifolds/L2Sphere/L2SphereVariable.cpp Manifolds/L2Sphere/L2SphereVector.cpp 
Manifolds/LowRank/LowRank.cpp Manifolds/LowRank/LowRankVariable.cpp Manifolds/LowRank/LowRankVector.cpp 
Manifolds/Oblique/Oblique.cpp Manifolds/Oblique/ObliqueVariable.cpp Manifolds/Oblique/ObliqueVector.cpp 
Manifolds/OrthGroup/OrthGroup.cpp Manifolds/OrthGroup/OrthGroupVariable.cpp Manifolds/OrthGroup/OrthGroupVector.cpp 
Manifolds/PreShapeCurves/PSCVariable.cpp Manifolds/PreShapeCurves/PSCVector.cpp Manifolds/PreShapeCurves/PreShapeCurves.cpp 
Manifolds/SPDManifold/SPDManifold.cpp Manifolds/SPDManifold/SPDVariable.cpp Manifolds/SPDManifold/SPDVector.cpp 
Manifolds/SPDTensor/SPDTVariable.cpp Manifolds/SPDTensor/SPDTVector.cpp Manifolds/SPDTensor/SPDTensor.cpp 
Manifolds/Sphere/Sphere.cpp Manifolds/Sphere/SphereVariable.cpp Manifolds/Sphere/SphereVector.cpp 
Manifolds/Stiefel/StieVariable.cpp Manifolds/Stiefel/StieVector.cpp Manifolds/Stiefel/Stiefel.cpp 
Others/ForDebug.cpp Others/MinPNormConHull.cpp Others/MyMatrix.cpp Others/Spline.cpp Others/Timer.cpp Others/randgen.cpp 
Others/SparseBLAS/nist_spblas.cpp 
Problems/Problem.cpp 
#Problems/juliaProblem.cpp 
#Problems/mexProblem.cpp 
Problems/ElasticCurvesRO/DriverElasticCurvesRO.cpp Problems/ElasticCurvesRO/ElasticCurvesRO.cpp 
Problems/EucFrechetMean/EucFrechetMean.cpp 
Problems/EucPosSpCd/EucPosSpCd.cpp 
Problems/EucQuadratic/EucQuadratic.cpp 
Problems/GrassRQ/GrassRQ.cpp 
Problems/KarcherMean/KarcherMean.cpp 
Problems/LRMatrixCompletion/LRMatrixCompletion.cpp 
Problems/ObliqueSparsePCA/ObliqueSparsePCA.cpp 
Problems/ObliqueTestSparsePCA/ObliqueTestSparsePCA.cpp 
Problems/OrthBoundingBox/OrthBoundingBox.cpp 
Problems/PreShapePathStraighten/PreShapePathStraighten.cpp 
Problems/SPDMean/SPDMean.cpp 
Problems/SPDTensorDL/SPDTensorDL.cpp 
Problems/ShapePathStraighten/ShapePathStraighten.cpp 
Problems/SphereConvexHull/SphereConvexHull.cpp 
Problems/StieBrockett/StieBrockett.cpp 
Problems/StieSoftICA/StieSoftICA.cpp 
Problems/StieSparseBrockett/StieSparseBrockett.cpp 
Problems/StieSumBrockett/StieSumBrockett.cpp 
Problems/WeightedLowrank/WeightedLowRank.cpp 
Solvers/LRBFGS.cpp Solvers/LRBFGSLPSub.cpp Solvers/LRTRSR1.cpp Solvers/MRankAdaptive.cpp Solvers/QuasiNewton.cpp Solvers/RBFGS.cpp Solvers/RBFGSLPSub.cpp Solvers/RBroydenFamily.cpp Solvers/RCG.cpp Solvers/RGS.cpp Solvers/RNewton.cpp Solvers/RSD.cpp Solvers/RTRNewton.cpp Solvers/RTRSD.cpp Solvers/RTRSR1.cpp Solvers/RWRBFGS.cpp Solvers/Solvers.cpp Solvers/SolversLS.cpp Solvers/SolversLSLPSub.cpp Solvers/SolversTR.cpp)
target_include_directories(roptlib PUBLIC ${ROPTLIB_INCLUDE_DIRECTORIES})
target_link_libraries(roptlib ${LAPACK} ${BLAS} ${M})
if(${CODE_PROFILING})
set_target_properties(roptlib PROPERTIES COMPILE_FLAGS "-pg -g" LINK_FLAGS "-pg -g")
endif()

# COMPILE EXECUTABLES

# Compile the simple test program from ROPTLIB
#add_executable(TestSimpleExample test/TestSimpleExample.cpp)
#target_link_libraries(TestSimpleExample roptlib)

# Compile the simple product test program from ROPTLIB
#add_executable(TestProduct test/TestProduct.h test/TestProduct.cpp)
#target_link_libraries(TestProduct roptlib)

# Compile the simple product test program from ROPTLIB
#add_executable(TestProductBrockett test/TestProductExample.cpp test/TestProductExample.cpp)
#target_link_libraries(TestProductBrockett roptlib)

export(TARGETS roptlib FILE ROPTLIB.cmake)
