
#include "test/TestMyMatrix.h"

using namespace ROPTLIB;

/*If the file is not compiled in Matlab and TESTSTIEBROCKETT is defined in def.h file, then using the following
main() function as the entrance. */
#if !defined(MATLAB_MEX_FILE) && defined(TESTMYMATRIX)

int main(void)
{
//	testEigenSymmetricM();
//	testExpSymmetricM();
	testLogSymmetricM();
	return 0;
}

void testEigenSymmetricM(void)
{
	double *M = new double[16 + 4 + 16];
	double *eigvalues = M + 16;
	double *eigvectors = eigvalues + 4;
	for (integer i = 0; i < 4; i++)
	{
		for (integer j = i; j < 4; j++)
		{
			M[i + j * 4] = genrandnormal();
			M[j + i * 4] = M[i + j * 4];
		}
	}
	Matrix A(M, 4, 4), E(eigvalues, 4, 1), V(eigvectors, 4, 4);
	Matrix::EigenSymmetricM(GLOBAL::U, A, E, V);
	delete[] M;
};

void testExpSymmetricM(void)
{
	double *M = new double[16 + 16];
	double *ExpM = M + 16;
	for (integer i = 0; i < 4; i++)
	{
		for (integer j = i; j < 4; j++)
		{
			M[i + j * 4] = genrandnormal();
			M[j + i * 4] = M[i + j * 4];
		}
	}
	Matrix A(M, 4, 4), B(ExpM, 4, 4);
	Matrix::ExpSymmetricM(GLOBAL::U, A, B);
	delete[] M;
};

void testLogSymmetricM(void)
{
	integer N = 4;
	double *M = new double[16 + 16 + 16];
	double *MMt = M + 16;
	double *LogM = MMt + 16;
	for (integer i = 0; i < 4; i++)
	{
		for (integer j = i; j < 4; j++)
		{
			M[i + j * 4] = genrandnormal();
			M[j + i * 4] = M[i + j * 4];
		}
	}
	dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, M, &N, M, &N, &GLOBAL::DZERO, MMt, &N);
	Matrix A(M, 4, 4), B(LogM, 4, 4), C(MMt, 4, 4);
	Matrix::LogSymmetricM(GLOBAL::U, C, B);
	delete[] M;
};

#endif
