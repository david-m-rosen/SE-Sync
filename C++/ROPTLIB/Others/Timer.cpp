
#include "Others/Timer.h"

/*Define the namespace*/
namespace ROPTLIB{

	unsigned long getTickCount(void)
	{
		//	return clock();

		unsigned long currentTime = 0;
//#ifdef MATLAB_MEX_FILE
//		mxArray *lhs[1], *rhs[1];
//		rhs[0] = mxCreateString("toc");
//		mexCallMATLAB(1, lhs, 1, rhs, "feval");
//		double result = mxGetScalar(lhs[0]);
//		return static_cast<unsigned long> (result * CLK_PS);
//#endif

#ifdef _WIN64
		LARGE_INTEGER li;
		QueryPerformanceFrequency(&li);
		long long dff = li.QuadPart;
		QueryPerformanceCounter(&li);
		currentTime = static_cast<unsigned long> (li.QuadPart * 1000000 / dff);
		//currentTime = GetTickCount() * 1000;
#elif _WIN32
		LARGE_INTEGER li;
		QueryPerformanceFrequency(&li);
		long long dff = li.QuadPart;
		QueryPerformanceCounter(&li);
		currentTime = static_cast<unsigned long> (li.QuadPart * 1000000 / dff);
		//currentTime = GetTickCount() * 1000;
#elif __APPLE__
#include "TargetConditionals.h"
#if TARGET_OS_MAC
		struct timeval current;
		gettimeofday(&current, NULL);
		currentTime = current.tv_sec * 1000000 + current.tv_usec;
#endif
#elif __linux
		struct timeval current;
		gettimeofday(&current, NULL);
		currentTime = current.tv_sec * 1000000 + current.tv_usec;
#endif // end of checking platforms

		return currentTime;
	}
}; /*end of ROPTLIB namespace*/
