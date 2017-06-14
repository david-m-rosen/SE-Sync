//
//  KMVariable.cpp
//  Updated
//
//  Created by Yaqing You on 5/4/16.
//  Copyright © 2016 Yaqing You. All rights reserved.
//

#include "Manifolds/ElasticShape/ShapeVariable.h"

/*Define the namespace*/
namespace ROPTLIB{
    ShapeVariable::ShapeVariable(integer r, integer c)
    {
        Element::Initialization(2, r, c);
    }

    ShapeVariable *ShapeVariable::ConstructEmpty(void) const
    {
        return new ShapeVariable(size[0], size[1]);
    }

}; /*end of ROPTLIB namespace*/
