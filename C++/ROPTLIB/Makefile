# Makefile for ROPTLIB. Test on ubuntu 16.04 LTS

# set compiler
CC = g++
CXXFLAGS:=-O3 -ffastmath -march=native -ggdb3

# default test problem is the Brockett cost function on the Stiefel manifold
TP?=TestSimpleExample

#the path of ROPTLIB
ROOTPATH:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))

# set the path of Julia
JULIA_DIR:=/home/whuang/Documents/julia

# directories of ROPTLIB header files
INCDIRS = -I$(ROOTPATH)/
INCDIRS += -I$(ROOTPATH)/BinaryFiles/
INCDIRS += -I$(ROOTPATH)/Julia/
INCDIRS += -I$(ROOTPATH)/Julia/useless/
INCDIRS += -I$(ROOTPATH)/Manifolds/
INCDIRS += -I$(ROOTPATH)/Manifolds/CpxNStQOrth/
INCDIRS += -I$(ROOTPATH)/Manifolds/ElasticShape/
INCDIRS += -I$(ROOTPATH)/Manifolds/EucPositive/
INCDIRS += -I$(ROOTPATH)/Manifolds/Euclidean/
INCDIRS += -I$(ROOTPATH)/Manifolds/Grassmann/
INCDIRS += -I$(ROOTPATH)/Manifolds/L2Sphere/
INCDIRS += -I$(ROOTPATH)/Manifolds/LowRank/
INCDIRS += -I$(ROOTPATH)/Manifolds/Oblique/
INCDIRS += -I$(ROOTPATH)/Manifolds/OrthGroup/
INCDIRS += -I$(ROOTPATH)/Manifolds/PreShapeCurves/
INCDIRS += -I$(ROOTPATH)/Manifolds/SPDManifold/
INCDIRS += -I$(ROOTPATH)/Manifolds/SPDTensor/
INCDIRS += -I$(ROOTPATH)/Manifolds/Sphere/
INCDIRS += -I$(ROOTPATH)/Manifolds/Stiefel/
INCDIRS += -I$(ROOTPATH)/Others/
INCDIRS += -I$(ROOTPATH)/Others/SparseBLAS/
INCDIRS += -I$(ROOTPATH)/Problems/
INCDIRS += -I$(ROOTPATH)/Problems/ElasticCurvesRO/
INCDIRS += -I$(ROOTPATH)/Problems/EucFrechetMean/
INCDIRS += -I$(ROOTPATH)/Problems/EucPosSpCd/
INCDIRS += -I$(ROOTPATH)/Problems/EucQuadratic/
INCDIRS += -I$(ROOTPATH)/Problems/GrassRQ/
INCDIRS += -I$(ROOTPATH)/Problems/KarcherMean/
INCDIRS += -I$(ROOTPATH)/Problems/LRMatrixCompletion/
INCDIRS += -I$(ROOTPATH)/Problems/ObliqueSparsePCA/
INCDIRS += -I$(ROOTPATH)/Problems/ObliqueTestSparsePCA/
INCDIRS += -I$(ROOTPATH)/Problems/OrthBoundingBox/
INCDIRS += -I$(ROOTPATH)/Problems/PreShapePathStraighten/
INCDIRS += -I$(ROOTPATH)/Problems/SPDMean/
INCDIRS += -I$(ROOTPATH)/Problems/SPDTensorDL/
INCDIRS += -I$(ROOTPATH)/Problems/ShapePathStraighten/
INCDIRS += -I$(ROOTPATH)/Problems/SphereConvexHull/
INCDIRS += -I$(ROOTPATH)/Problems/StieBrockett/
INCDIRS += -I$(ROOTPATH)/Problems/StieSoftICA/
INCDIRS += -I$(ROOTPATH)/Problems/StieSparseBrockett/
INCDIRS += -I$(ROOTPATH)/Problems/StieSumBrockett/
INCDIRS += -I$(ROOTPATH)/Problems/WeightedLowrank/
INCDIRS += -I$(ROOTPATH)/Solvers/
INCDIRS += -I$(ROOTPATH)/cwrapper/
INCDIRS += -I$(ROOTPATH)/cwrapper/blas/
INCDIRS += -I$(ROOTPATH)/cwrapper/lapack/
INCDIRS += -I$(ROOTPATH)/test/
# ROPTLIB C++ files
CPPS += $(ROOTPATH)/Manifolds/Element.cpp $(ROOTPATH)/Manifolds/LinearOPE.cpp $(ROOTPATH)/Manifolds/Manifold.cpp $(ROOTPATH)/Manifolds/ProductElement.cpp $(ROOTPATH)/Manifolds/ProductManifold.cpp $(ROOTPATH)/Manifolds/SharedSpace.cpp $(ROOTPATH)/Manifolds/SmartSpace.cpp 
CPPS += $(ROOTPATH)/Manifolds/CpxNStQOrth/CSOVariable.cpp $(ROOTPATH)/Manifolds/CpxNStQOrth/CSOVector.cpp $(ROOTPATH)/Manifolds/CpxNStQOrth/CpxNStQOrth.cpp 
CPPS += $(ROOTPATH)/Manifolds/ElasticShape/ElasticShape.cpp $(ROOTPATH)/Manifolds/ElasticShape/ShapeVariable.cpp $(ROOTPATH)/Manifolds/ElasticShape/ShapeVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/EucPositive/EucPosVariable.cpp $(ROOTPATH)/Manifolds/EucPositive/EucPosVector.cpp $(ROOTPATH)/Manifolds/EucPositive/EucPositive.cpp 
CPPS += $(ROOTPATH)/Manifolds/Euclidean/EucVariable.cpp $(ROOTPATH)/Manifolds/Euclidean/EucVector.cpp $(ROOTPATH)/Manifolds/Euclidean/Euclidean.cpp 
CPPS += $(ROOTPATH)/Manifolds/Grassmann/GrassVariable.cpp $(ROOTPATH)/Manifolds/Grassmann/GrassVector.cpp $(ROOTPATH)/Manifolds/Grassmann/Grassmann.cpp 
CPPS += $(ROOTPATH)/Manifolds/L2Sphere/L2Sphere.cpp $(ROOTPATH)/Manifolds/L2Sphere/L2SphereVariable.cpp $(ROOTPATH)/Manifolds/L2Sphere/L2SphereVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/LowRank/LowRank.cpp $(ROOTPATH)/Manifolds/LowRank/LowRankVariable.cpp $(ROOTPATH)/Manifolds/LowRank/LowRankVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/Oblique/Oblique.cpp $(ROOTPATH)/Manifolds/Oblique/ObliqueVariable.cpp $(ROOTPATH)/Manifolds/Oblique/ObliqueVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/OrthGroup/OrthGroup.cpp $(ROOTPATH)/Manifolds/OrthGroup/OrthGroupVariable.cpp $(ROOTPATH)/Manifolds/OrthGroup/OrthGroupVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/PreShapeCurves/PSCVariable.cpp $(ROOTPATH)/Manifolds/PreShapeCurves/PSCVector.cpp $(ROOTPATH)/Manifolds/PreShapeCurves/PreShapeCurves.cpp 
CPPS += $(ROOTPATH)/Manifolds/SPDManifold/SPDManifold.cpp $(ROOTPATH)/Manifolds/SPDManifold/SPDVariable.cpp $(ROOTPATH)/Manifolds/SPDManifold/SPDVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/SPDTensor/SPDTVariable.cpp $(ROOTPATH)/Manifolds/SPDTensor/SPDTVector.cpp $(ROOTPATH)/Manifolds/SPDTensor/SPDTensor.cpp 
CPPS += $(ROOTPATH)/Manifolds/Sphere/Sphere.cpp $(ROOTPATH)/Manifolds/Sphere/SphereVariable.cpp $(ROOTPATH)/Manifolds/Sphere/SphereVector.cpp 
CPPS += $(ROOTPATH)/Manifolds/Stiefel/StieVariable.cpp $(ROOTPATH)/Manifolds/Stiefel/StieVector.cpp $(ROOTPATH)/Manifolds/Stiefel/Stiefel.cpp 
CPPS += $(ROOTPATH)/Others/ForDebug.cpp $(ROOTPATH)/Others/MinPNormConHull.cpp $(ROOTPATH)/Others/MyMatrix.cpp $(ROOTPATH)/Others/Spline.cpp $(ROOTPATH)/Others/Timer.cpp $(ROOTPATH)/Others/randgen.cpp 
CPPS += $(ROOTPATH)/Others/SparseBLAS/nist_spblas.cpp 
CPPS += $(ROOTPATH)/Problems/Problem.cpp $(ROOTPATH)/Problems/juliaProblem.cpp $(ROOTPATH)/Problems/mexProblem.cpp 
CPPS += $(ROOTPATH)/Problems/ElasticCurvesRO/DriverElasticCurvesRO.cpp $(ROOTPATH)/Problems/ElasticCurvesRO/ElasticCurvesRO.cpp 
CPPS += $(ROOTPATH)/Problems/EucFrechetMean/EucFrechetMean.cpp 
CPPS += $(ROOTPATH)/Problems/EucPosSpCd/EucPosSpCd.cpp 
CPPS += $(ROOTPATH)/Problems/EucQuadratic/EucQuadratic.cpp 
CPPS += $(ROOTPATH)/Problems/GrassRQ/GrassRQ.cpp 
CPPS += $(ROOTPATH)/Problems/KarcherMean/KarcherMean.cpp 
CPPS += $(ROOTPATH)/Problems/LRMatrixCompletion/LRMatrixCompletion.cpp 
CPPS += $(ROOTPATH)/Problems/ObliqueSparsePCA/ObliqueSparsePCA.cpp 
CPPS += $(ROOTPATH)/Problems/ObliqueTestSparsePCA/ObliqueTestSparsePCA.cpp 
CPPS += $(ROOTPATH)/Problems/OrthBoundingBox/OrthBoundingBox.cpp 
CPPS += $(ROOTPATH)/Problems/PreShapePathStraighten/PreShapePathStraighten.cpp 
CPPS += $(ROOTPATH)/Problems/SPDMean/SPDMean.cpp 
CPPS += $(ROOTPATH)/Problems/SPDTensorDL/SPDTensorDL.cpp 
CPPS += $(ROOTPATH)/Problems/ShapePathStraighten/ShapePathStraighten.cpp 
CPPS += $(ROOTPATH)/Problems/SphereConvexHull/SphereConvexHull.cpp 
CPPS += $(ROOTPATH)/Problems/StieBrockett/StieBrockett.cpp 
CPPS += $(ROOTPATH)/Problems/StieSoftICA/StieSoftICA.cpp 
CPPS += $(ROOTPATH)/Problems/StieSparseBrockett/StieSparseBrockett.cpp 
CPPS += $(ROOTPATH)/Problems/StieSumBrockett/StieSumBrockett.cpp 
CPPS += $(ROOTPATH)/Problems/WeightedLowrank/WeightedLowRank.cpp 
CPPS += $(ROOTPATH)/Solvers/LRBFGS.cpp $(ROOTPATH)/Solvers/LRBFGSLPSub.cpp $(ROOTPATH)/Solvers/LRTRSR1.cpp $(ROOTPATH)/Solvers/MRankAdaptive.cpp $(ROOTPATH)/Solvers/QuasiNewton.cpp $(ROOTPATH)/Solvers/RBFGS.cpp $(ROOTPATH)/Solvers/RBFGSLPSub.cpp $(ROOTPATH)/Solvers/RBroydenFamily.cpp $(ROOTPATH)/Solvers/RCG.cpp $(ROOTPATH)/Solvers/RGS.cpp $(ROOTPATH)/Solvers/RNewton.cpp $(ROOTPATH)/Solvers/RSD.cpp $(ROOTPATH)/Solvers/RTRNewton.cpp $(ROOTPATH)/Solvers/RTRSD.cpp $(ROOTPATH)/Solvers/RTRSR1.cpp $(ROOTPATH)/Solvers/RWRBFGS.cpp $(ROOTPATH)/Solvers/Solvers.cpp $(ROOTPATH)/Solvers/SolversLS.cpp $(ROOTPATH)/Solvers/SolversLSLPSub.cpp $(ROOTPATH)/Solvers/SolversTR.cpp 

# convert a string to upper case.
UPPER_TP  = $(shell echo $(TP) | tr a-z A-Z)

# make a binary file, which is called in command line
ROPTLIB:
	$(CC) -O3 -w -std=c++0x $(ROOTPATH)/test/$(TP).cpp $(CPPS) $(INCDIRS) -D$(UPPER_TP) -llapack -lblas -lm -o $(TP)

#make a library
libropt.so:
	$(CC) -w -std=c++0x -shared -fPIC -O3 $(MANIFOLDS) $(OTHERS) $(PROBLEMS) $(SOLVERS) $(INCDIRS) -llapack -lblas -lm -o $@

JULIA_LIB:=$(JULIA_DIR)/usr/lib
JULIA_SRC:=$(JULIA_DIR)/src
JULIA_INC:=$(JULIA_DIR)/usr/include
CPPFLAGS:=-I$(JULIA_INC) -I$(JULIA_SRC) -I$(JULIA_SRC)/support
LDFLAGS:=-L$(JULIA_LIB)
LDLIBS=-ljulia
export LD_LIBRARY_PATH:=$(JULIA_LIB):$(JULIA_LIB)/julia

# make a shared library, which is used by Julia
JuliaROPTLIB:
	$(CC) -O3 -shared -fPIC -std=c++0x $(ROOTPATH)/test/$(TP).cpp $(CPPS) $(INCDIRS) -D$(UPPER_TP) $(CPPFLAGS) $(LDFLAGS) -Wl,-rpath,$(JULIA_LIB) -lm $(LDLIBS) -DJULIA_LIB_DIR=\"$(JULIA_DIR)/lib/julia\" -llapack -lblas -o $(TP).so
