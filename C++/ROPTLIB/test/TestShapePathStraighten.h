//
//  TestShapePathStraighten.hpp
//  Updated
//
//  Created by Yaqing You on 2/8/16.
//  Copyright © 2016 Yaqing You. All rights reserved.
//

#ifndef TestShapePathStraighten_h
#define TestShapePathStraighten_h

//#include <stdio.h>

#include <fstream>
//#include "Others/ForDebug.h"  //cannot include fordebug here (cause the "abs" problems)
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
#include "Problems/ShapePathStraighten/ShapePathStraighten.h"

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
#include "Problems/mexProblem.h"   //

#include "Others/def.h"

#include "Problems/ElasticCurvesRO/DriverElasticCurvesRO.h"

using namespace ROPTLIB;

void testShapePathStraighten();


#endif /* TestShapePathStraighten_h */
