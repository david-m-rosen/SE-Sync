/*
compute the minimum p-norm in a convex hull:
min_{x \in conv{v1, v2, ... v_p}} \|x\|_P,
where $P$ is given by a function handle, v_i, i=1,...,p are tangent vectors.

If pure C++ is used, then the algorithm stated in [SY1992] or Riemannian method can used.
If Matlab is used, then the matlab function "quadprog" is used.

[SY1992]: K. Sekitami and Y. Yamamoto, A recursive algorithm for finding the minimum
norm point in a polytope and a pair of closest points in two polytopes.

-----WH
*/

#ifndef MINPNORMCONHULL_H
#define MINPNORMCONHULL_H

#include "Problems/Problem.h"
#include "Manifolds/Manifold.h"
#include "Solvers/QuasiNewton.h"
#include "Problems/SphereConvexHull/SphereConvexHull.h"
#include "Others/MyMatrix.h"
#include "Others/def.h"

#define RIEMANNIANCONHULL // RECURSIVEMETHOD // 

/*Define the namespace*/
namespace ROPTLIB{

	/*Compute min_{y in convex hull of gfs, and gfs are tangent vectors at the tangent space at x} ||y||_2 */
	extern double MinPNormConHull(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, Vector *Soln, double *YtY, integer inc);

	/*Compute min_{y in convex hull of gfs, and gfs are tangent vectors at the tangent space at x} ||y|| */
	extern double MinPNormConHull(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, QuasiNewton *solver, 
		void (QuasiNewton::*Hv)(Vector *v, Vector *result), Vector *Soln, double *YtPY, integer inc);

	/*Compute min_{y in convex hull of gfs, and gfs are tangent vectors at the tangent space at x} ||y||_P */
	extern double MinPNormConHullRecursive(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, 
		QuasiNewton *solver, void (QuasiNewton::*Hv)(Vector *v, Vector *result), Vector *Soln, double *YtPY, integer inc);

	extern double MinPNormConHullRecursivesubprob(const Manifold *Mani, Variable *x, Vector **Ys, Vector **PYs, integer LYs, bool *idxYs, double *Pnorm2Ys, double NumErr, Vector *initialX,
		QuasiNewton *solver, void (QuasiNewton::*Hv)(Vector *v, Vector *result), Vector *Soln, integer &times);

	/*Compute min_{y in convex hull of gfs, and gfs are tangent vectors at the tangent space at x} ||y|| */
	extern double MinPNormConHullRMethod(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, QuasiNewton *solver, 
		void (QuasiNewton::*Hv)(Vector *v, Vector *result), Vector *Soln, double *YtPY, integer inc);

	/*Compute min_{y in convex hull of gfs, and gfs are tangent vectors at the tangent space at x} ||y|| */
	extern double MinPNormConHullMatlab(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, QuasiNewton *solver, 
		void (QuasiNewton::*Hv)(Vector *v, Vector *result), Vector *Soln, double *YtPY, integer inc);

}; /*end of ROPTLIB namespace*/

#endif
