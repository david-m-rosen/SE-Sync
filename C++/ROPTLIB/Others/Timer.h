#ifndef TIMER_H
#define TIMER_H

#include <ctime>

#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif

#ifdef _WIN64
#undef VOID
#include <windows.h>
#elif _WIN32
#undef VOID
#include <windows.h>
#elif __APPLE__
    #include "TargetConditionals.h"
    #if TARGET_OS_MAC
    #include <unistd.h>
    #include <sys/time.h>
    #include <netinet/in.h>
    #endif
#elif __linux
#include <unistd.h>
#include <sys/time.h>
#include <netinet/in.h>
#endif // end of checking platforms

//#include "def.h"

#define CLK_PS 1000000
//#define CLK_PS CLOCKS_PER_SEC

/*Define the namespace*/
namespace ROPTLIB{

	unsigned long getTickCount();
}; /*end of ROPTLIB namespace*/
#endif
