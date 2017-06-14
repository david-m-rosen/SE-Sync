//
//  KMVector.hpp
//  Updated
//
//  Created by Yaqing You on 5/4/16.
//  Copyright © 2016 Yaqing You. All rights reserved.
//

#ifndef SHAPEVECTOR_H
#define SHAPEVECTOR_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{
    class ShapeVector : public Element{
    public:
        /* r=numP, c=dim */
        ShapeVector(integer r, integer c);
        
        virtual ShapeVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif /* SHAPEVECTOR_H */


