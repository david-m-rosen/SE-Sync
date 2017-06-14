
#include "test/TestCSO.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTCSO)

std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testCSO();
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

#endif

void testCSO()
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
	genrandseed(tt);
	integer n = 6, p = 2;
	// Obtain an initial iterate by taking the Q factor of qr decomposition
	// lapack is used
	CSOVariable SCOX(n, p);

	//double *ptr = SCOX.ObtainWriteEntireData();
	//double a = 1.0 / sqrt(8) * 100000;
	//for (integer i = 0; i < 8; i++)
	//	ptr[i] = a;
	SCOX.RandInManifold();

	// Define the manifold
	CpxNStQOrth Domain(n, p);
	//Domain.SetHasHHR(true);

	Domain.CheckParams();

	Domain.CheckIntrExtr(&SCOX);
	Domain.CheckRetraction(&SCOX);
	Domain.CheckDiffRetraction(&SCOX);
	Domain.CheckIsometryofVectorTransport(&SCOX);
	Domain.CheckIsometryofInvVectorTransport(&SCOX);
	Domain.CheckLockingCondition(&SCOX);
	Domain.CheckcoTangentVector(&SCOX);
	Domain.CheckVecTranComposeInverseVecTran(&SCOX);
	Domain.CheckTranHInvTran(&SCOX);
	Domain.CheckHaddScaledRank1OPE(&SCOX);
}

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	testCSO();
}

#endif
