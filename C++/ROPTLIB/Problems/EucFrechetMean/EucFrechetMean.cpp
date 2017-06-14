
#include "Problems/EucFrechetMean/EucFrechetMean.h"

/*Define the namespace*/
namespace ROPTLIB{

	EucFrechetMean::EucFrechetMean(double *W, double *D, integer num, integer dim)
	{
		Weights = W;
		Data = D;
		Dim = dim;
		Num = num;
	};

	EucFrechetMean::~EucFrechetMean(void)
	{
	};

	double EucFrechetMean::f(Variable *x) const
	{
		EucVariable *Eucx = dynamic_cast<EucVariable *> (x);
		double result = 0;
		const double *xM = Eucx->ObtainReadData();
		for (integer i = 0; i < Num; i++) // LAPACK
		{
			for (integer j = 0; j < Dim; j++)
				result += Weights[i] * (xM[j] - Data[j + i * Dim]) * (xM[j] - Data[j + i * Dim]);
		}
		return result;
	};

	void EucFrechetMean::Grad(Variable *x, Vector *gf) const
	{
		EucVariable *Eucx = dynamic_cast<EucVariable *> (x);
		EucVector *Eucgf = dynamic_cast<EucVector *> (gf);
		const double *xM = Eucx->ObtainReadData();
		double *gfTV = Eucgf->ObtainWriteEntireData();

		for (integer i = 0; i < Dim; i++) // LAPACK
		{
			gfTV[i] = 0;
			for (integer j = 0; j < Num; j++)
			{
				gfTV[i] += 2.0 * Weights[j] * (xM[i] - Data[i + j * Dim]);
			}
		}
	};

	void EucFrechetMean::HessianEta(Variable *x, Vector *etax, Vector *xix) const
	{
		double sum = 0;
		for (integer i = 0; i < Num; i++)
		{
			sum += 2 * Weights[i];
		}
		Domain->ScaleTimesVector(x, sum, etax, xix);
	};
}; /*end of ROPTLIB namespace*/
