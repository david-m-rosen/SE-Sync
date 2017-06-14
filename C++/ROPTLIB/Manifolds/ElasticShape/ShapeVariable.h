//
//  KMVariable.hpp
//  Updated
//
//  Created by Yaqing You on 5/4/16.
//  Copyright © 2016 Yaqing You. All rights reserved.
//

#ifndef ShapeVariable_H
#define ShapeVariable_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{
    class ShapeVariable : public Element{
    public:
        // r : number of points on a curve
        // l : dimension of a space the curves are in
        
        /* Construct an empty variable on the space with only size information */
        ShapeVariable(integer r, integer l);
        
        /* Generate an initial curve */
        //void Generate();
        
        /* Creat an object of KMVariable with same size as this KMVariable */
        virtual ShapeVariable *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif /* ShapeVariable_H */
