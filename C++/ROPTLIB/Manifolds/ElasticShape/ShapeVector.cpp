//
//  KMVector.cpp
//  Updated
//
//  Created by Yaqing You on 5/4/16.
//  Copyright © 2016 Yaqing You. All rights reserved.
//

#include "Manifolds/ElasticShape/ShapeVector.h"

/*Define the namespace*/
namespace ROPTLIB{
    ShapeVector::ShapeVector(integer r, integer c)
    {
        Element::Initialization(2, r, c);
    }

    ShapeVector *ShapeVector::ConstructEmpty(void) const
    {
        return new ShapeVector(size[0], size[1]);
	}
}; /*end of ROPTLIB namespace*/
