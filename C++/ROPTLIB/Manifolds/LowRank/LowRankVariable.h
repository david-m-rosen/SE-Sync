/*
This file defines the class of a point on the low-rank manifold R_r^{m times n}, which is represented by 
Gr(r, m) times R^{r times r} times Gr(r, n).

SmartSpace --> ProductElement --> LowRankVariable

---- WH
*/

#ifndef LOWRANKVARIABLE_H
#define LOWRANKVARIABLE_H

#include "Manifolds/ProductElement.h"
#include "Manifolds/Grassmann/GrassVariable.h"
#include "Manifolds/Euclidean/EucVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LowRankVariable : public ProductElement{
	public:
		/*Construct an empty variable on the low-rank manifold R_r^{m times n} with only size information.*/
		LowRankVariable(integer m, integer n, integer r);

		/*Destruct by deleting all variables*/
		virtual ~LowRankVariable(void);

		/*Create an object of LowRankVariable with same size as this LowRankVariable.*/
		virtual LowRankVariable *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of LOWRANKVARIABLE_H
