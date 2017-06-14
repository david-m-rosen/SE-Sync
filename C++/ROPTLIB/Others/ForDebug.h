/*
This file defines some functions to help debugging.
In the current version, it only contains on function which prints a double array.

-----WH
*/

#ifndef FORDEBUG_H
#define FORDEBUG_H

#include "Others/MyMatrix.h"
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class ForDebug{
	public:
		static void Print(const char *name, const double *M, integer row, integer col = 1, integer num = 1);
		static double NormF(const double *V, integer length);
	};
}; /*end of ROPTLIB namespace*/
#endif // end of FORDEBUG_H
