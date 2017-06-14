/*
This file defines the class for all the product of manifolds
It defines the common properties and features of all the product of manifolds.
Users can write their own product of manifolds by deriving from this class.
Note that if a function requires some arguments can be the same,
then users need to guarantee that the derived function also
support this property.

Manifold --> ProductManifold

---- WH
*/

#ifndef PRODUCTMANIFOLD_H
#define PRODUCTMANIFOLD_H

/*This is used to help to check whether there is a memory leakage or not.*/
#define CHECKMEMORY

/*ProductVariable and ProductVector are just ProductElement*/
#define ProdVariable ProductElement
#define ProdVector ProductElement

#include "Manifolds/Manifold.h"
#include "Manifolds/ProductElement.h"
#include "Manifolds/Euclidean/EucVariable.h"
#include "Manifolds/Euclidean/EucVector.h"
#include <cstdarg>
#include <map>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class ProductManifold : public Manifold{
	public:
		/*Constructor of ProductManifold.
		An example of using this constructor to generate St(2,3)^2 \times Euc(2) is:
		integer n = 3, p = 2, m = 2;
		integer numofmanis = 2;
		integer numofmani1 = 2;
		integer numofmani2 = 1;
		Stiefel mani1(n, p);
		Euclidean mani2(m);
		ProductManifold ProdMani(numofmanis, &mani1, numofmani1, &mani2, numofmani2);
		See examples in TestProduct.cpp and more in test files.	*/
		ProductManifold(integer numberofmanifolds, ...);

		/*Constructor of ProductManifold
			An example of using this constructor to generate St(2,3)^2 \times Euc(2) is:
			integer n = 3, p = 2, m = 2;
			Stiefel mani1(n, p);
			Euclidean mani2(m);
			Manifold **inmanifolds = new Manifold* [2]; inmanifolds[0] = &mani1; inmanifolds[1] = &mani2;
			integer inpowsinterval = {0, 2, 3};
			ProductManifold ProdMani(inmanifolds, 2, inpowsinterval, 3);
			* The first argument indicates that there are two kinds of manifolds St(2, 3) and Euc(2).
			* The second argement indicates that the length of "inmanifolds" is 2.
			* The third argument indicates the number for each manifold. The number of St(2, 3) is inpowsinterval[1] - inpowsinterval[0] = 2
			and the number of Euc(2) is inpowsinterval[2] - inpowsinterval[1] = 1. Therefore, the product manifold if St(2, 3)^2 \times Euc(2).
			* The fourth argument indicates the number of all manifolds, i.e., 2 Stiefel + 1 Euclidean = 3. */
		ProductManifold(Manifold **inmanifolds, integer innumofmani, integer *inpowsinterval, integer innumoftotalmani);

		/*In ProductManifold, true of IsIntrApproach means that the representation of each component of tangent vector of
		each manifold is specified by IsIntrApproach parameter of each manifolds. True of IsIntrApproach does NOT mean
		the representation of each component of tangent vector uses intrinsic representation. */
		void SetEMPTYINTR();

		/*Destructor of ProductManifold. Delete EMPTYINTR, EMPTYEXTR, manifolds, and powsinterval;*/
		virtual ~ProductManifold();

		/*Define the Riemannian metric: g_x(etax, xix).
		The default one is summation of the metrics of all manifolds, i.e.,
		g_x(etax, xix) = Sum g_{x_i}(etax_i, xix_i), where x_i, etax_i, and xix_i are components of x, etax, and xix respectively. */
		virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;

		/*Define the action of the linear operator, result = Hx(etax), where Hx is an linear operator on T_x M and
		etax is a tangent vector in T_x M. Note that etax and result must be different, i.e., calling this
		function by
		LinearOPEEta(x, Hx, result, result);
		is illegal.
		The default one is the matrix multiplication, i.e., result = Hx * etax*/
		virtual void LinearOPEEta(Variable *x, LinearOPE *Hx, Vector *etax, Vector *result) const;

		/*Compute result = scalar1 * etax + scalar2 * xix; etax and result can be a same argument or xix and result can be a same argument, i.e.,
		calling this function by
		VectorLinearCombination(x, scalar1, etax, scalar2, result, result); or VectorLinearCombination(x, scalar1, result, scalar2, xix, result);
		is legal.
		However, VectorLinearCombination(x, scalar1, result, scalar2, result, result) is illegal.*/
		virtual void VectorLinearCombination(Variable *x, double scalar1, Vector *etax, double scalar2, Vector *xix, Vector *result) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax
		The component of v and results use the representation specified by each manifold if IsIntrApproach is true.
		Otherwise, the component of v and results use the extrinsic representation.
		etax and result can be a same argument, i.e.,
		calling this function by
		Projection(x, result, result);
		is legal.
		Default function: Let x = (x_1, dots, x_n), v=(v_1, dots, v_n), result = (r_1, dots, r_n).
		The default computation is r_i = P_{T_{x_i} M_i} v_i for all i. */
		virtual void Projection(Variable *x, Vector *v, Vector *result) const;

		/*Randomly generate N vectors in the tangent space of x. The resulting vectors are stored in result_arr.
		This is used for the gradient sampling method, which is not supported yet.*/
		virtual void RandomTangentVectors(Variable *x, integer N, Vector **result_arr) const;

		/*Compute the retraction result = R_x(etax). A stepsize is also input for information.
		Default: Let x = (x_1, dots, x_n), etax=(etax_1, dots, etax_n), result = (r_1, dots, r_n).
		r_i = R_{x_i}(eta_i), for all i. */
		virtual void Retraction(Variable *x, Vector *etax, Variable *result, double instepsize) const;

		/*Compute the retraction result = R_x(etax).
		Default: Let x = (x_1, dots, x_n), etax=(etax_1, dots, etax_n), result = (r_1, dots, r_n).
		r_i = R_{x_i}(eta_i), for all i. */
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

		/*Compute the tangent vector result satisfying
		g_y(\mathcal{T}_{R_etax}(xix), xiy) = g_x(xix, result) for all xix \in T_x M,
		where y = R_x(etax), xiy \in T_y M.
		This cotangent vector is used in the RBFGS defined in [RW2012]
		[RW2012]: W. Ring and B. Wirth. Optimization methods on Riemannian manifolds and their application to shape space.
		SIAM Journal on Optimization, 22(2):596?27, January 2012
		xiy and result can be a same argument, i.e.,
		calling this function by
		coTangentVector(x, etax, y, result, result);
		is legal.
		Default: compute the cotangent vector for each component */
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Computes the vector transport by differentiated retraction, i.e., result = \mathcal{T}_{R_etax} (xix)
		if the input IsEtaXiSameDir is true, then this means etax and xix are along a same direction. This implies
		the computations may be simplied.
		xix and result can be a same argument, i.e.,
		calling this function by
		DiffRetraction(x, etax, y, result, result);
		is legal.
		Default: compute the vector transport by differentiated retraction for each component */
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*computes beta = \|etax\| / \|\mathcal{T}_{R_etax} etax\|
		Default: etax = (etax_1, dots, eta_n);
		beta = \sqrt(sum \|eta_i\|_2^2) / \sqrt(sum \|\mathcal{T}_{R_{etax_i}} etax_i\|) */
		virtual double Beta(Variable *x, Vector *etax) const;

		/*compute the distance between two points on the manifold.
		Default: sqrt(dist_M1(x1_1, x2_2)^2 + dist_M2(x1_2, x2_2)^2 + ... + dist_M1(x1_n, x2_n)^2)
		where x1_i is the i-th component of the product element.*/
		virtual double Dist(Variable *x1, Variable *x2) const;

		/*Computes the vector transport, i.e., result = \mathcal{T}_etax (xix)
		xix and result can be a same argument, i.e.,
		calling this function by
		VectorTransport(x, etax, y, result, result);
		is legal.
		Default: compute the vector transport for each component */
		virtual void VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

		/*Computes the inverse vector transport, i.e., result = \mathcal{T}_etax^{-1} (xiy)
		xiy and result can be a same argument, i.e.,
		calling this function by
		InverseVectorTransport(x, etax, y, result, result);
		is legal.
		Default: compute the inverse vector transport for each component */
		virtual void InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Compute result = diag(\mathcal{T}_1, dots, \mathcal{T}_n) * H * diag(\mathcal{T}_1^{-1}, dots, \mathcal{T}_n^{-1}).
		Hx and result can be a same argument, i.e.,
		calling this function by
		TranHInvTran(x, etax, y, result, result);
		is legal.
		Default: This is implemented by call the functions "HInvTran" and "TranH" in class "Manifold". */
		virtual void TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;

		/*Compute result = Hx + scalar * etax * xix^{\flat}.
		Hx and result can be a same argument, i.e.,
		calling this function by
		HaddScaledRank1OPE(x, result, scalar, etax, xix, result);
		is legal.
		Default: compute xix^{\flat} = (xix_1^\flat, cdots, xix_n^\flat) first, then compute esult = Hx + scalar * etax * xix^{\flat}*/
		virtual void HaddScaledRank1OPE(Variable *x, LinearOPE *Hx, double scalar, Vector *etax, Vector *xix, LinearOPE *result) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax;
		For this function, the components of v and result are represented by extrinsic representations.
		etax and result can be a same argument, i.e.,
		calling this function by
		ExtrProjection(x, result, result);
		is legal.
		Default function: Let x = (x_1, dots, x_n), v=(v_1, dots, v_n), result = (r_1, dots, r_n).
		The default computation is r_i = P_{T_{x_i} M_i} v_i for all i. */
		virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;


		/*Each component of etax must use extrinsic representation.
		This function turn the extrinsic representation of a component of etax to intrinsic representation if the
		component manifold has true IsIntrApproach. Otherwise, the component of etax still use extrinsic representation.
		See details in [Hua2013, Section 9.5]
		[Hua2013]: W. Huang. Optimization algorithms on Riemannian manifolds with applications.
		PhD thesis, Florida State University, Department of Mathematics, 2013
		etax and result must be different arguments, i.e.,
		calling this function by
		ObtainIntr(x, result, result);
		is illegal.
		Default: If the "IsIntrApproach" in M_i is true, then call M_i::ObtainIntr(x_i, etax_i, result_i).
		Otherwise, result_i <-- etax_i.*/
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*Each component of etax use representation specified by the "IsIntrApproach" of each manifold.
		This function turn the intrinsic representation of a component of etax to extrinsic representation if the
		component manifold has true IsIntrApproach.
		See details in [Hua2013, Section 9.5]
		[Hua2013]: W. Huang. Optimization algorithms on Riemannian manifolds with applications.
		PhD thesis, Florida State University, Department of Mathematics, 2013
		etax and result must be different arguments, i.e.,
		calling this function by
		ObtainExtr(x, result, result);
		is illegal.
		Default: If the "IsIntrApproach" in M_i is true, then call M_i::ObtainExtr(x_i, intretax_i, result_i).
		Otherwise, result_i <-- etax_i.*/
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

		/*Compute the Riemannian gradient from the Euclidean gradient of a function;
		The function is defined in "prob".
		egf is the Euclidean gradient; result is the output which is the Riemannian gradient.
		egf and result can be a same argument, i.e.,
		calling this function by
		EucGradToGrad(x, result, result, prob);
		is legal.
		Default: call the EucGradToGrad for each component manifold. */
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*Compute the Riemannian action of Hessian from the Euclidean action of Hessian of a function;
		The function is defined in "prob".
		the input exix is exix = EucHess [etax], where EucHess is the Euclidean Hessian,
		the output result is result = RieHess[eta], where RieHess is the Riemannian Hessian.
		exix and result can be a same argument, i.e.,
		calling this function by
		EucHvToHv(x, etax, result, result, prob);
		is legal.
		Default: call the EucHvToHv for each component manifold. */
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*Get the idx-th kind of manifolds.*/
		inline Manifold *GetManifold(integer idx) const { if (idx < numofmani) return manifolds[idx]; return nullptr; };

	protected:
		/*Turn ProdElement to Element. It is not the same as dynamic_cast. This function has not been used yet.
		I don't remember why I wrote it*/
		void ProdElementToElement(const ProductElement *ProdElem, Element *Elem) const;

		/*Turn ProdElement to Element. It is not the same as dynamic_cast. This function has not been used yet.
		I don't remember why I wrote it*/
		void ElementToProdElement(const Element *Elem, ProductElement *ProdElem) const;

		Manifold **manifolds; /*Store all kinds of manifolds*/
		integer numofmani; /*The number of kinds of manifolds, i.e., the length of manifolds*/
		integer *powsinterval; /*Each manifold are located in what interval*/
		integer numoftotalmani; /*The number of all manifolds.*/

		/*An example for store St(2, 3) ^ 2 \times Euc(2) is given below:
		(Not exactly C++ code. just give an idea)
		Stiefel mani1(3, 2);
		Euclidean mani2(2);
		manifolds = {&mani1, &mani2}; // two kinds of manifolds
		numofmani = 2; // two kinds of manifold
		powsinterval = {0, 2, 3}; // from 0 to 2-1 is St(2,3); from 2 to 3-1 is Euc(2);
		numoftotalmani = 3; // 2 Stiefel + 1 Euclidean = 3, i.e., 3 manifolds in total.
		*/
	};
}; /*end of ROPTLIB namespace*/
#endif // end of PRODUCTMANIFOLD_H