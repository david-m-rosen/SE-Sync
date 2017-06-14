//
//  TestPreShapePathStraighten.hpp
//  headertry
//
//  Created by Yaqing You on 11/3/15.
//  Copyright © 2015 Yaqing You. All rights reserved.
//

#ifndef TestPreShapePathStraighten_h
#define TestPreShapePathStraighten_h

#include <fstream>

#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>
#include "test/DriverMexProb.h"

#include "Manifolds/PreShapeCurves/PSCVariable.h"
#include "Manifolds/PreShapeCurves/PSCVector.h"


#include "Manifolds/PreShapeCurves/PreShapeCurves.h"
#include "Problems/PreShapePathStraighten/PreShapePathStraighten.h"

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

#include "Problems/ElasticCurvesRO/DriverElasticCurvesRO.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTPRESHAPEPATHSTRAIGHTEN)
int main(void);
#endif

void testPreShapePathStraighten();


#endif /* TestPreShapePathStraighten_h */
