
#include "Problems/EucQuadratic/EucQuadratic.h"

/*Define the namespace*/
namespace ROPTLIB{

	EucQuadratic::EucQuadratic(double *M, integer dim)
	{
		Dim = dim;
		A = M;
	};

	EucQuadratic::~EucQuadratic(void)
	{
	};

	double EucQuadratic::f(Variable *x) const
	{
		const double *v = x->ObtainReadData();
		SharedSpace *Temp = new SharedSpace(1, Dim);
		double *temp = Temp->ObtainWriteEntireData();

		char *transn = const_cast<char *> ("n");
		double one = 1, zero = 0;
		integer inc = 1, N = Dim;
		// temp <- A * v, details: http://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
		dgemv_(transn, &N, &N, &one, A, &N, const_cast<double *> (v), &inc, &zero, temp, &inc);

		x->AddToTempData("Ax", Temp);

		// output v^T temp, details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
		return ddot_(&N, const_cast<double *> (v), &inc, temp, &inc);
	};

	void EucQuadratic::Grad(Variable *x, Vector *gf) const
	{
		double *gfTV = gf->ObtainWriteEntireData();
		const SharedSpace *Temp = x->ObtainReadTempData("Ax");
		const double *v = Temp->ObtainReadData();

		integer N = Dim, inc = 1;
		// gfTV <- v, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&N, const_cast<double *> (v), &inc, gfTV, &inc);
		double two = 2;
		// gfTV <- 2 * gfTV, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dscal_(&N, &two, gfTV, &inc);
	};

	void EucQuadratic::HessianEta(Variable *x, Vector *etax, Vector *xix) const
	{
		const double *v = etax->ObtainReadData();
		double *xixTV = xix->ObtainWriteEntireData();

		char *transn = const_cast<char *> ("n");
		integer N = Dim, inc = 1;
		double two = 2, zero = 0;
		// xixTV <- 2 * A * v, details: http://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
		dgemv_(transn, &N, &N, &two, A, &N, const_cast<double *> (v), &inc, &zero, xixTV, &inc);
	};
}; /*end of ROPTLIB namespace*/
