/*
This file defines the class for temporary storage. An object of SharedSpace type can be 
attached to a point on a manifold or a tangent vector. This is usually used to avoid
redundent computations.

SmartSpace --> SharedSpace

---- WH
*/

#ifndef SHAREDSPACE_H
#define SHAREDSPACE_H

#include "Manifolds/SmartSpace.h"
#include "Manifolds/Element.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*Declaration of Element. Element has been defined somewhere.
	Element is the based class of points in manifolds, tangent vectors in tangent spaces.
	The SharedSpace can be attached to an object of Element. */
	class Element;

	class SharedSpace : public SmartSpace{
	public:

		/*There are two kinds of constructors for SharedSpace.
		(1), SharedSpace(integer numberofdimensions, ...);
		(2), SharedSpace(Element *inelement);
		The first one uses an array of double numbers to store temporary data.
		The second one uses an object of Element to store temporary data.
		*/

		/*Initialize the SharedSpace. If one want to create a 10 by 3 by 2 by 5 tensor to store
		temporary data, then call this function by
		SharedSpace temp(4, 10, 3, 2, 5);
		The first argument indicates the number of following parameters.
		The next 4 arguments means the dimensions of the tensor.
		This SharedSpace will not allocate memory in this function. */
		SharedSpace(integer numberofdimensions, ...);

		/*Initialize the SharedSpace. If one want to store a gradient, which is an object of class derived from Element,
		then call this function by
		SharedSpace temp(Gradient)*/
		SharedSpace(Element *inelement);

		/*Destruct*/
		~SharedSpace();

		/*Create an object of SharedSpace with same size as this ShareSpace. The space is not allocated at this stage.*/
		virtual SharedSpace *ConstructEmpty(void) const;

		/*Copy this SharedSpace to "eta" SharedSpace. After calling this function,
		this SharedSpace and "eta" SharedSpace will use same space to store data. */
		virtual void CopyTo(SharedSpace *eta) const;

		/*Print the data. The string "name" is to mark the output such that user can find the output easily.*/
		virtual void Print(const char *name = "") const;

		/*Get the shared Element. It is useful only the temporary data is an object of Element, not an array of double numbers.*/
		Element *GetSharedElement(void) const;
	private:
		/*The shared Element*/
		Element *SharedElement;
	};
}; /*end of ROPTLIB namespace*/

#endif // end of SHAREDSPACE_H
