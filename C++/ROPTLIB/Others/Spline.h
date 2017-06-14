/*
This file defines function for computing cubic spline curve for various functions.
TODO: Define the class to be a spline curve, not just a class contains many static functions.

-----WH
*/

#ifndef SPLINE_H
#define SPLINE_H

//#include <cmath>

#include <iostream>
#include <limits>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class Spline{
	public:
		static int SplineUniformPeriodic(const double *Y, int n, double h, double *coefs); // periodic boundary condition: for closed curves
		static int SplinePeriodic(const double *X, const double *Y, int n, double *coefs);

		static int SplineUniformSlopes(const double *Y, int n, double h, double *coefs); // Hermite boundary conditions: for open curves (derivative-free from)
		static int SplineSlopes(const double *X, const double *Y, int n, double *coefs);

		static int SolveTridiagonalSystem(double *d, double *ud, double *ld, double *vec, double *s, int n);
		static int SolvePeriodicSystem(double *d, double *ud, double *ld, double *vec, double *s, int nn);

		static double ValSpline(const double *coefs, const double *breaks, int N, double t);
		static double ValSplineUniform(const double *coefs, int N, double h, double t);
		static void FirstDeri(const double *coefs, int N, double *dericoefs);
		static double ValFirstDeriUniform(const double *dericoefs, int N, double h, double t);
		static double ValFirstDeri(const double *dericoefs, const double *breaks, int N, double t);
		static void SecondDeri(const double *coefs, int N, double *dericoefs);
		static double ValSecondDeriUniform(const double *dericoefs, int N, double h, double t);
		static double ValSecondDeri(const double *dericoefs, const double *breaks, int N, double t);
	};
}; /*end of ROPTLIB namespace*/
#endif // SPLINE_H
