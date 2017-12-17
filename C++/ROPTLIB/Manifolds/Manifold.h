/*
This file defines the abstract base class for all the manifolds
It defines the common properties and features of all the manifolds.
Users can write their own manifolds by deriving from this class.
Note that if a function requires some arguments can be the same,
then users need to guarantee that the derived function also
support this property.

Manifold

---- WH
*/

#ifndef MANIFOLD_H
#define MANIFOLD_H


//#include <cmath>

#include "Problems/Problem.h"
#include "Manifolds/Element.h"
#include "Manifolds/LinearOPE.h"
#include "Manifolds/SharedSpace.h"
#include <iomanip>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*declaration of Problem and ProductManifold. Manifold class should know
	the classes Problem and ProductManifold have been defined somewhere.*/
	class Problem;
	class ProductManifold;

	class Manifold{
	public:
		/*Indicate this class is abstract*/
		virtual ~Manifold(void) = 0;

		/*Define the Riemannian metric: g_x(etax, xix).
		The default one is the Euclidean metric.*/
		virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;

		/*Define the action of the linear operator, result = Hx(etax), where Hx is an linear operator on T_x M and
		etax is a tangent vector in T_x M. Note that etax and result must be different, i.e., calling this
		function by
		LinearOPEEta(x, Hx, result, result);
		is illegal.
		The default one is the matrix multiplication, i.e., result = Hx * etax*/
		virtual void LinearOPEEta(Variable *x, LinearOPE *Hx, Vector *etax, Vector *result) const;

		/*Compute result = scalar * etax; etax and result can be a same argument, i.e.,
		calling this function by
		ScaleTimesVector(x, scalar, result, result);
		is legal. */
		virtual void ScaleTimesVector(Variable *x, double scalar, Vector *etax, Vector *result) const;

		/*Compute result = etax + xix; etax and result can be a same argument or xix and result can be a same argument, i.e.,
		calling this function by
		VectorAddVector(x, etax, result, result); or VectorAddVector(x, result, xix, result);
		is legal.
		However, VectorAddVector(x, result, result, result) is illegal.*/
		virtual void VectorAddVector(Variable *x, Vector *etax, Vector *xix, Vector *result) const;

		/*Compute result = etax - xix; etax and result can be a same argument or xix and result can be a same argument, i.e.,
		calling this function by
		VectorMinusVector(x, etax, result, result); or VectorMinusVector(x, result, xix, result);
		is legal.
		However, VectorMinusVector(x, result, result, result) is illegal.*/
		virtual void VectorMinusVector(Variable *x, Vector *etax, Vector *xix, Vector *result) const;

		/*Compute result = scalar * etax + xix; etax and result can be a same argument or xix and result can be a same argument, i.e.,
		calling this function by
		scalarVectorAddVector(x, scalar, etax, result, result); or scalarVectorAddVector(x, scalar, result, xix, result);
		is legal.
		However, scalarVectorAddVector(x, scalar, result, result, result) is illegal.*/
		virtual void scalarVectorAddVector(Variable *x, double scalar, Vector *etax, Vector *xix, Vector *result) const;

		/*Compute result = scalar * etax - xix; etax and result can be a same argument or xix and result can be a same argument, i.e.,
		calling this function by
		scalarVectorMinusVector(x, scalar, etax, result, result); or scalarVectorMinusVector(x, scalar, result, xix, result);
		is legal.
		However, scalarVectorMinusVector(x, scalar, result, result, result) is illegal.*/
		virtual void scalarVectorMinusVector(Variable *x, double scalar, Vector *etax, Vector *xix, Vector *result) const;

		/*Compute result = scalar1 * etax + scalar2 * xix; etax and result can be a same argument or xix and result can be a same argument, i.e.,
		calling this function by
		VectorLinearCombination(x, scalar1, etax, scalar2, result, result); or VectorLinearCombination(x, scalar1, result, scalar2, xix, result);
		is legal.
		However, VectorLinearCombination(x, scalar1, result, scalar2, result, result) is illegal.*/
		virtual void VectorLinearCombination(Variable *x, double scalar1, Vector *etax, double scalar2, Vector *xix, Vector *result) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax
		etax and result can be a same argument, i.e.,
		calling this function by
		Projection(x, result, result);
		is legal.
		Default: result <-- etax*/
		virtual void Projection(Variable *x, Vector *etax, Vector *result) const;

		/*Randomly generate N vectors in the tangent space of x. The resulting vectors are stored in result_arr.
		This is used for the gradient sampling method, which is not supported yet.*/
		virtual void RandomTangentVectors(Variable *x, integer N, Vector **result_arr) const;

		/*Compute the retraction result = R_x(etax). A stepsize is also input for information.
		Default: result = x + etax;*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result, double instepsize) const;

		/*Compute the retraction result = R_x(etax).
		Default: result = x + etax;*/
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
			Default: result <-- xiy */
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Computes the vector transport by differentiated retraction, i.e., result = \mathcal{T}_{R_etax} (xix)
		if the input IsEtaXiSameDir is true, then this means etax and xix are along a same direction. This implies
		the computations may be simplied.
		xix and result can be a same argument, i.e.,
		calling this function by
		DiffRetraction(x, etax, y, result, result);
		is legal.
		Default: result <-- xix */
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*computes beta = \|etax\| / \|\mathcal{T}_{R_etax} etax\|
		Default: beta <-- 1*/
		virtual double Beta(Variable *x, Vector *etax) const;

		/*compute the distance between two points on the manifold.
		Default: \|x1 - x2\|_F*/
		virtual double Dist(Variable *x1, Variable *x2) const;

		/*Computes the vector transport, i.e., result = \mathcal{T}_etax (xix)
		xix and result can be a same argument, i.e.,
		calling this function by
		VectorTransport(x, etax, y, result, result);
		is legal.
		Default: result <-- xix*/
		virtual void VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

		/*Computes the inverse vector transport, i.e., result = \mathcal{T}_etax^{-1} (xiy)
		xiy and result can be a same argument, i.e.,
		calling this function by
		InverseVectorTransport(x, etax, y, result, result);
		is legal.
		Default: result <-- xiy*/
		virtual void InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Compute result = H(:, start : end) * \mathcal{T}^{-1}, where H(:, start : end) denotes the matrix formed by columns from "start" to "end".
		The movitation to use "start" and "end" is for the product manifolds. See function TranHInvTran in ProductManifold.h.
		Hx and result can be a same argument, i.e.,
		calling this function by
		HInvTran(x, etax, y, result, start, end, result);
		is legal.
		Default: result <-- Hx*/
		virtual void HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*Compute result = \mathcal{T} * H(start : end, :), where H(start : end, :) denotes the matrix formed by cows from "start" to "end".
		The movitation to use "start" and "end" is for the product manifolds. See function TranHInvTran in ProductManifold.h.
		Hx and result can be a same argument, i.e.,
		calling this function by
		TranH(x, etax, y, result, start, end, result);
		is legal.
		Default: result <-- Hx*/
		virtual void TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*Compute result = \mathcal{T} * H * \mathcal{T}^{-1}.
		Hx and result can be a same argument, i.e.,
		calling this function by
		TranHInvTran(x, etax, y, result, result);
		is legal.
		Default: result <-- Hx*/
		virtual void TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;

		/*Compute result = Hx + scalar * etax * xix^{\flat}.
		Hx and result can be a same argument, i.e.,
		calling this function by
		HaddScaledRank1OPE(x, result, scalar, etax, xix, result);
		is legal.
		Default: result = Hx + scaler * etax * xix^T*/
		virtual void HaddScaledRank1OPE(Variable *x, LinearOPE *Hx, double scalar, Vector *etax, Vector *xix, LinearOPE *result) const;

		/*Compute etaxflat = etax^{\flat}.
		etax and etaxflat can be a same argument, i.e.,
		calling this function by
		ObtainEtaxFlat(x, etaxflat, etaxflat);
		is legal.
		Default: etaxflat <-- etax*/
		virtual void ObtainEtaxFlat(Variable *x, Vector *etax, Vector *etaxflat) const;

		/*Compute the intrinsic representation of a tangent vector etax,
		See details in [Hua2013, Section 9.5]
		[Hua2013]: W. Huang. Optimization algorithms on Riemannian manifolds with applications.
		PhD thesis, Florida State University, Department of Mathematics, 2013
		etax and result must be different arguments, i.e.,
		calling this function by
		ObtainIntr(x, result, result);
		is illegal.
		Default: result <-- etax*/
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*Compute the extrinsic representation of a tangent vector etax,
		See details in [Hua2013, Section 9.5]
		[Hua2013]: W. Huang. Optimization algorithms on Riemannian manifolds with applications.
		PhD thesis, Florida State University, Department of Mathematics, 2013
		etax and result must be different arguments, i.e.,
		calling this function by
		ObtainIntr(x, result, result);
		is illegal.
		Default: result <-- intretax */
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;


		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax;
		For this function, both etax and result are represented by intrinsic representations.
		etax and result can be a same argument, i.e.,
		calling this function by
		IntrProjection(x, result, result);
		is legal.
		Default: result <-- etax*/
		virtual void IntrProjection(Variable *x, Vector *etax, Vector *result) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax;
		For this function, both etax and result are represented by extrinsic representations.
		etax and result can be a same argument, i.e.,
		calling this function by
		ExtrProjection(x, result, result);
		is legal.
		Default: result <-- etax*/
		virtual void ExtrProjection(Variable *x, Vector *etax, Vector *result) const;

		/*Get the name of current manifold*/
		inline std::string GetName(void) const { return name; };

		/*Get the dimension of the manifold*/
		inline integer GetIntrDim(void) const { return IntrinsicDim; };

		/*Get the dimension of the ambient space*/
		inline integer GetExtrDim(void) const { return ExtrinsicDim; };

		/*Set the representation of tangent vectors. True if intrinsic is used, false otherwise*/
		inline void SetIsIntrApproach(bool IsIntrAppr) const { IsIntrApproach = IsIntrAppr; };

		/*Get the representation of tangent vectors. True if intrinsic is used, false otherwise*/
		inline bool GetIsIntrinsic(void) const { return IsIntrApproach; };

		/*Get an empty tangent vector using intrinisic representation*/
		inline const Vector *GetEMPTYINTR(void) const { return EMPTYINTR; };

		/*Get an empty taangent vector using extrinisic representation*/
		inline const Vector *GetEMPTYEXTR(void) const { return EMPTYEXTR; };

		/*Check whether the locking condition is satisfied or not.*/
		inline bool GetHasLockCon(void) const { return HasLockCon; };

		/*Set whether users want to apply the idea in [HGA2015, Section 4.1] to the vector transport defined in the function "VectorTransport".
			[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
			SIAM Journal on Optimization, 25(3):1660?685,2015.	.*/
		inline void SetHasHHR(bool inHasHHR) { HasHHR = inHasHHR; };

		/*Check whether the idea in [HGA2015, Section 4.1] is used or not*/
		inline bool GetHasHHR(void) { return HasHHR; };

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*Compute the Riemannian gradient from the Euclidean gradient of a function;
		The function is defined in "prob".
		egf is the Euclidean gradient; result is the output which is the Riemannian gradient.
		egf and result can be a same argument, i.e.,
		calling this function by
		EucGradToGrad(x, result, result, prob);
		is legal.
		It is a pure virtual function. It must be overloaded by derived class */
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *result, const Problem *prob) const = 0;

		/*Compute the Riemannian action of Hessian from the Euclidean action of Hessian of a function;
		The function is defined in "prob".
		the input exix is exix = EucHess [etax], where EucHess is the Euclidean Hessian,
		the output result is result = RieHess[eta], where RieHess is the Riemannian Hessian.
		exix and result can be a same argument, i.e.,
		calling this function by
		EucHvToHv(x, etax, result, result, prob);
		is legal.
		It is a pure virtual function.It must be overloaded by derived class */
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* result, const Problem *prob) const = 0;

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/* ===============
		The functions below are used to check the correctness of the Riemannian functions.
		They can be used for people who need to write their own manifold.
		================= */

		/*Check the correctness of the functions: "ObtainIntr" and "ObtainExtr",
		under the assumption that function "ExtrProjection" is correct.
		If the metrics for extrinsic representation and intrinsic representation have been implemented,
		then the outputs of "extr inp" equivalent to "intr inp" means the vector transport is isometric.*/
		virtual void CheckIntrExtr(Variable *x) const;

		/*Check the correctness of the Retraction function: "Retraction",
		under the assumption that functions "ExtrProjection", "ObtainIntr", and "ObtainExtr" are correct.*/
		virtual void CheckRetraction(Variable *x) const;

		/*Check the correctness of the Retraction function: "DiffRetraction",
		under the assumption that functions "ExtrProjection", "ObtainIntr", "ObtainExtr" and "Retraction" are correct.*/
		virtual void CheckDiffRetraction(Variable *x, bool IsEtaXiSameDir = true) const;

		/*Check whether the vector transport and differentiated retraction satisfy the locking condition,
			where the locking condition is define in [HGA2015, (2.8)]
			[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
			SIAM Journal on Optimization, 25(3):1660?685, 2015
			under the assumption that functions "ExtrProjection", "ObtainIntr", "ObtainExtr", "Retraction" and "DiffRetraction" are correct.*/
		virtual void CheckLockingCondition(Variable *x) const;

		/*Check the correctness of the cotangent vector function: "coTangentVector",
		under the assumption that functions "ExtrProjection", "ObtainIntr", "ObtainExtr", "Retraction" and "DiffRetraction" are correct.*/
		virtual void CheckcoTangentVector(Variable *x) const;

		/*Check the isometry of vector transport. This function gives correct result only when extrinsic representation is used.
		If intrinsic representation is used, please use "CheckIntrExtr" instead. */
		virtual void CheckIsometryofVectorTransport(Variable *x) const;

		/*Check the isometry of inverse vector transport. This function gives correct result only when extrinsic representation is used.
		If intrinsic representation is used, please use "CheckIntrExtr" instead. */
		virtual void CheckIsometryofInvVectorTransport(Variable *x) const;

		/*Check whether \mathcal{T} \circ \mathcal{T}^{-1} is identity or not, where \mathcal{T} is defined by function "VectorTransport"
		and \mathcal{T}^{-1} is defined by function "InverseVectorTransport".*/
		virtual void CheckVecTranComposeInverseVecTran(Variable *x) const;

		/*Check whether TranHInvTran is correct or not.*/
		virtual void CheckTranHInvTran(Variable *x) const;

		/*Check whether "HaddScaledRank1OPE" is correct or not*/
		virtual void CheckHaddScaledRank1OPE(Variable *x) const;

	protected:
		mutable bool HasHHR;/*Mark whether the idea in [HGA2015, Section 4.1] is used or not*/
		mutable bool IsIntrApproach; /*Mark whether intrinsic representation is used for tangent vector or not*/
		bool UpdBetaAlone; /*For some manifolds, the vector transport and retraction already satisfy the locking condition except a correct beta.
						   By setting this parameter true, the algorithm does not use the Householer transformations in [HGA2015, Section 4.1].
						   Only computing the beta can make the locking condition is satisfied.*/
		std::string name; /*The name of this manifold*/
		integer IntrinsicDim; /*The dimension of this manifold*/
		integer ExtrinsicDim; /*The dimension of the ambient space.*/
		Vector *EMPTYINTR; /*The empty tangent vector using intrinsic representation*/
		Vector *EMPTYEXTR; /*The empty tangent vector using extrinsic representation*/

		bool HasLockCon;/*Mark whether the locking condition is satisfied or not.*/

		/* The next 5 functions modified above 5 functions such that the locking condition is satisfied.
		The idea in [HGA2015, Section 4.1] is used.
		They are not required to be satisfied.*/

		/*Apply idea in [HGA2015, Section 4.1] to the function "VectorTransport". */
		virtual void LCVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

		/*Apply idea in [HGA2015, Section 4.1] to the function "LCInverseVectorTransport". */
		virtual void LCInverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Apply idea in [HGA2015, Section 4.1] to the function "LCHInvTran". */
		virtual void LCHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*Apply idea in [HGA2015, Section 4.1] to the function "LCTranH". */
		virtual void LCTranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*Apply idea in [HGA2015, Section 4.1] to the function "LCTranHInvTran". */
		virtual void LCTranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;

		/*The function computes unit vectors used in above LC* functions. The idea is in [HGA2015, Section 4.1]. */
		virtual void Obtainnu1nu2forLC(Variable *x, Vector *etax, Variable *y) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of MANIFOLD_H
