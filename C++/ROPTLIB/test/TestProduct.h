/*
This is the test file to check the correctenss of classes: ProductElement and ProductManifold.

---- WH
*/

#ifndef TESTPRODUCT_H
#define TESTPRODUCT_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>

#include "Manifolds/Euclidean/EucVariable.h"
#include "Manifolds/Euclidean/EucVector.h"
#include "Problems/EucFrechetMean/EucFrechetMean.h"
#include "Problems/EucQuadratic/EucQuadratic.h"

#include "Problems/StieBrockett/StieBrockett.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/Stiefel.h"

#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"

#include "Solvers/SolversTR.h"
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

#include "Others/def.h"

#include "Manifolds/ProductElement.h"
#include "Manifolds/ProductManifold.h"

#include "Manifolds/L2Sphere/L2Sphere.h"
#include "Manifolds/L2Sphere/L2SphereVariable.h"
#include "Manifolds/L2Sphere/L2SphereVector.h"

#include "Manifolds/Oblique/Oblique.h"
#include "Manifolds/Oblique/ObliqueVariable.h"
#include "Manifolds/Oblique/ObliqueVector.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTPRODUCT)
int main(void);
#endif

void testProduct(void);

#endif // end of TESTPRODUCT_H
