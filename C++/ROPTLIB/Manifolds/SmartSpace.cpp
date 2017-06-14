
#include "Manifolds/SmartSpace.h"

/*Define the namespace*/
namespace ROPTLIB{

	void SmartSpace::Initialization(integer numberofdimensions, ...)
	{
		va_list argptr;
		//	va_start(argptr, numberofdimensions);
		ls = numberofdimensions;
		size = new integer[ls];
		va_start(argptr, numberofdimensions);
		for (integer i = 0; i < ls; i++)
			size[i] = va_arg(argptr, integer);
		va_end(argptr);

		length = 1;
		for (integer i = 0; i < ls; i++)
			length *= size[i];
		Space = nullptr;
		sharedtimes = nullptr;
	};

	void SmartSpace::NewMemoryOnWrite(void)
	{
		if (sharedtimes == nullptr)
		{
			NewMemory();
			sharedtimes = new integer;
			*sharedtimes = 1;
		}
		else
			if ((*sharedtimes) > 1)
			{
				NewMemory();
				(*sharedtimes)--;
				sharedtimes = new integer;
				*sharedtimes = 1;
			}
	};

	void SmartSpace::CopyOnWrite(void)
	{
		if (sharedtimes == nullptr)
		{
			NewMemory();
			sharedtimes = new integer;
			*sharedtimes = 1;
		}
		else
			if ((*sharedtimes) > 1)
			{
				double *ptr = Space;
				NewMemory();
				(*sharedtimes)--;
				sharedtimes = new integer;
				*sharedtimes = 1;

				integer N = length, inc = 1;
				// Space <- ptr, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
				dcopy_(&N, ptr, &inc, Space, &inc);
			}
	};

	void SmartSpace::RandUnform(double start, double end)
	{
		NewMemoryOnWrite();
		double ell = end - start;
		for (integer i = 0; i < length; i++)
			Space[i] = genrandreal() * ell + start;
	};

	void SmartSpace::RandGaussian(double mean, double variance)
	{
		NewMemoryOnWrite();
		for (integer i = 0; i < length; i++)
			Space[i] = (genrandnormal() + mean) * variance;
	};

	void SmartSpace::CopyTo(SmartSpace *eta) const
	{
		if (this == eta || eta->Space == Space)
			return;
		bool IsSameSize = true;

		if (eta->ls != ls)
		{
			IsSameSize = false;
		}
		else
		{
			for (integer i = 0; i < ls; i++)
			{
				if (eta->size[i] != size[i])
				{
					IsSameSize = false;
					break;
				}
			}
		}
		if (eta->sharedtimes != nullptr && *(eta->sharedtimes) == 1 && IsSameSize)
		{
			integer N = length, inc = 1;
			// (eta->Space) <- Space, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
			dcopy_(&N, Space, &inc, eta->Space, &inc);
			return;
		}

		if (eta->sharedtimes != nullptr && *(eta->sharedtimes) > 1)
		{
			(*(eta->sharedtimes))--;
		}
		else
			if (eta->sharedtimes != nullptr && *(eta->sharedtimes) == 1)
			{

#ifdef CHECKMEMORYDELETED
				(*CheckMemoryDeleted)[eta->sharedtimes] = *(eta->sharedtimes);
#endif
				delete eta->sharedtimes; eta->sharedtimes = nullptr;
				delete[] eta->Space; eta->Space = nullptr;
			}

		if (sharedtimes != nullptr)
			(*sharedtimes)++;
		eta->sharedtimes = sharedtimes;
		eta->Space = Space;

		if (eta->ls != ls)
		{
			delete[] eta->size;
			eta->size = new integer[ls];
			eta->ls = ls;
		}
		for (integer i = 0; i < ls; i++)
			eta->size[i] = size[i];
		eta->length = length;
	};

	void SmartSpace::NewMemory(void)
	{
		try{
			Space = new double[length];
		}
		catch (std::bad_alloc exception)
		{
			printf("Catch exception:%s\n", exception.what());
		}
	};

	const double *SmartSpace::ObtainReadData(void) const
	{
		return Space;
	};

	double *SmartSpace::ObtainWriteEntireData(void)
	{
		NewMemoryOnWrite();
		return Space;
	};

	double *SmartSpace::ObtainWritePartialData(void)
	{
		CopyOnWrite();
		return Space;
	};

	SmartSpace::~SmartSpace(void)
	{
		if (size != nullptr)
			delete[] size;

		if (sharedtimes != nullptr)
		{
#ifdef CHECKMEMORYDELETED
			(*CheckMemoryDeleted)[sharedtimes] = *sharedtimes;
#endif

			(*sharedtimes)--;
			if (*sharedtimes == 0 && Space != nullptr)
			{
				delete sharedtimes; sharedtimes = nullptr;
				delete[] Space; Space = nullptr;
			}
		}
	};

	void SmartSpace::Print(const char * name) const
	{
		integer product = 1;
		for (integer i = 2; i < ls; i++)
			product *= size[i];
		if (Space == nullptr)
		{
			if (size == nullptr)
			{
				printf("%s is an empty data with size 0", name);
			}
			else
			{
				printf("%s is an empty data with size %d", name, size[0]);
			}
			for (integer i = 1; i < ls; i++)
				printf(" x %d", size[i]);
			printf("\n");
		}
		else
			if (ls == 1 || (ls > 1 && size[1] * product == 1))
			{
				printf("%s , shared times:%d, shared times address:%p\n", name, *sharedtimes, sharedtimes);
				for (integer i = 0; i < length; i++)
					printf("%.10e\n", Space[i]);
			}
			else
				if (ls == 2 || product == 1)
				{
					printf("%s , shared times:%d, shared times address:%p\n", name, *sharedtimes, sharedtimes);
					for (integer j = 0; j < size[0]; j++)
					{
						for (integer k = 0; k < size[1]; k++)
						{
							printf("%.10e\t", Space[j + size[0] * k]);
						}
						printf("\n");
					}
				}
				else
				{
					integer st = 0, ed = 0, row = size[0], col = size[1];
					integer *idices = new integer[ls + 1];
					double *ptr = Space;
					for (integer i = 2; i < ls + 1; i++)
						idices[i] = 0;
					while (1)
					{
						printf("%s(:,:", name);
						for (integer i = 2; i < ls; i++)
							printf(",%d", idices[i]);
						printf("), shared times:%d\n", *sharedtimes);
						for (integer j = 0; j < row; j++)
						{
							for (integer k = 0; k < col; k++)
							{
								printf("%.10e\t", ptr[j + row * k]);
							}
							printf("\n");
						}
						ptr += row * col;
						idices[2]++;
						for (integer i = 2; i < ls; i++)
						{
							if (idices[i] == size[i])
							{
								idices[i] = 0;
								idices[i + 1]++;
							}
						}
						if (idices[ls] == 1)
							break;
					}
					delete[] idices;
				}
	};

	void SmartSpace::SetByParams(integer *insize, integer inls, integer inlength, integer *insharedtimes, double *inSpace)
	{
		size = insize;
		ls = inls;
		length = inlength;
		sharedtimes = insharedtimes;
		Space = inSpace;
	};

	void SmartSpace::DeleteBySettingNull(void)
	{
		size = nullptr;
		sharedtimes = nullptr;
		Space = nullptr;
	};
}; /*end of ROPTLIB namespace*/
