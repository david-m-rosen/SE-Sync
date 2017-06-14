
#include "Manifolds/SharedSpace.h"

/*Define the namespace*/
namespace ROPTLIB{

	SharedSpace::SharedSpace(integer numberofdimensions, ...)
	{
		ls = 1;
		va_list argptr;
		va_start(argptr, numberofdimensions);
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
		SharedElement = nullptr;
	};

	SharedSpace::SharedSpace(Element *inelement)
	{
		SharedElement = inelement;
		size = nullptr;
		ls = 0;
		length = 0;
		sharedtimes = nullptr;
		Space = nullptr;
	};

	SharedSpace::~SharedSpace(void)
	{
		if (SharedElement != nullptr)
			delete SharedElement;
	};

	SharedSpace *SharedSpace::ConstructEmpty(void) const
	{
		SharedSpace *ptr = new SharedSpace(1, 1);
		this->CopyTo(ptr);
		ptr->Space = nullptr;
		ptr->sharedtimes = nullptr;
		if (sharedtimes != nullptr)
		{
			(*sharedtimes)--;
		}
		return ptr;
	};

	void SharedSpace::CopyTo(SharedSpace *eta) const
	{
		SmartSpace::CopyTo(eta);
		if (SharedElement == nullptr && eta->SharedElement != nullptr)
		{
			delete eta->SharedElement;
			eta->SharedElement = nullptr;
		}
		if (SharedElement != nullptr && eta->SharedElement == nullptr)
		{
			eta->SharedElement = SharedElement->ConstructEmpty();
			SharedElement->CopyTo(eta->SharedElement);
		}
		if (SharedElement != nullptr && eta->SharedElement != nullptr)
		{
			SharedElement->CopyTo(eta->SharedElement);
		}
	};

	void SharedSpace::Print(const char *name) const
	{
		if (Space != nullptr)
		{
			SmartSpace::Print(name);
		}
		if (SharedElement != nullptr)
		{
			SharedElement->Print(name);
		}
	};

	Element *SharedSpace::GetSharedElement(void) const
	{
		return SharedElement;
	};
}; /*end of ROPTLIB namespace*/
