
#include "Manifolds/Element.h"

/*Define the namespace*/
namespace ROPTLIB{

	Element::~Element(void)
	{
		RemoveAllFromTempData();
	};

	void Element::CopyTo(Element *eta) const
	{
		SmartSpace::CopyTo(eta);
		MAP::const_iterator thisiter = TempData.begin();
		MAP::const_iterator etaiter, etaiterpre;
		for (thisiter = TempData.begin(); thisiter != TempData.end(); thisiter++)
		{
			etaiter = eta->TempData.find(thisiter->first);
			if (etaiter != eta->TempData.end())
			{
				thisiter->second->CopyTo(etaiter->second);
			}
			else
			{
				SharedSpace *Temp = thisiter->second->ConstructEmpty();
				thisiter->second->CopyTo(Temp);
				eta->AddToTempData(thisiter->first, Temp);
			}
		}
		if (TempData.size() < eta->TempData.size())
		{
			for (etaiter = eta->TempData.begin(); etaiter != eta->TempData.end();)
			{
				thisiter = TempData.find(etaiter->first);
				if (thisiter == TempData.end())
				{
					etaiterpre = etaiter;
					etaiter++;
					eta->RemoveFromTempData(etaiterpre->first);
				}
			}
		}
	};

	void Element::RandUnform(double start, double end)
	{
		RemoveAllFromTempData();
		SmartSpace::RandUnform(start, end);
	};

	void Element::RandGaussian(double mean, double variance)
	{
		RemoveAllFromTempData();
		SmartSpace::RandGaussian(mean, variance);
	};

	double *Element::ObtainWriteEntireData(void)
	{
		RemoveAllFromTempData();
		return SmartSpace::ObtainWriteEntireData();
	};

	double *Element::ObtainWritePartialData(void)
	{
		RemoveAllFromTempData();
		return SmartSpace::ObtainWritePartialData();
	};

	//void Element::NewMemoryOnWrite(void)
	//{
	//	RemoveAllFromTempData();
	//	SmartSpace::NewMemoryOnWrite();
	//};

	//void Element::CopyOnWrite(void)
	//{
	//	RemoveAllFromTempData();
	//	SmartSpace::CopyOnWrite();
	//};

	void Element::Print(const char *name, bool isonlymain) const
	{
		if (TempData.size() > 0 && !isonlymain)
			printf("=================Main data: %s=========================\n", name);
		SmartSpace::Print(name);

		if (TempData.size() > 0 && !isonlymain)
		{
			MAP::const_iterator thisiter;
			for (thisiter = TempData.begin(); thisiter != TempData.end(); thisiter++)
			{
				printf("=================Temp data in %s ================\n", name);
				thisiter->second->Print(thisiter->first.c_str());
			}
			printf("=================end of output: %s=========================\n", name);
		}
	};

	void Element::RandInManifold(void)
	{
		printf("Warning: RandInManifold has not been overloaded!\n");
	};

	void Element::AddToTempData(std::string name, SharedSpace * &Temp)
	{
		MAP::iterator thisiter;
		thisiter = TempData.find(name);
		if (thisiter == TempData.end())
		{
			TempData.insert(std::pair<std::string, SharedSpace *>(name, Temp));
		}
		else
		{
			Temp->CopyTo(thisiter->second);
			delete Temp;
		}
		Temp = nullptr;
	};

	const SharedSpace *Element::ObtainReadTempData(std::string name) const
	{
		MAP::const_iterator thisiter;
		thisiter = TempData.find(name);
		if (thisiter == TempData.end())
		{
			printf("Error: TempData %s does not exist!\n", name.c_str());
			return nullptr;
		}

		return thisiter->second;
	};

	SharedSpace *Element::ObtainWriteTempData(std::string name)
	{
		MAP::iterator thisiter;
		thisiter = TempData.find(name);
		if (thisiter == TempData.end())
		{
			printf("Error: TempData %s does not exist!\n", name.c_str());
			return nullptr;
		}

		return thisiter->second;
	};

	void Element::RemoveFromTempData(std::string name)
	{
		MAP::iterator thisiter;
		thisiter = TempData.find(name);
		if (thisiter != TempData.end())
		{
			delete thisiter->second;
			TempData.erase(thisiter);
		}
	};

	void Element::RemoveAllFromTempData(void)
	{
		MAP::iterator thisiter;
		for (thisiter = TempData.begin(); thisiter != TempData.end(); thisiter++)
		{
			delete thisiter->second;
		}
		TempData.clear();
	};

	bool Element::TempDataExist(std::string name) const
	{
		MAP::const_iterator thisiter;
		thisiter = TempData.find(name);
		if (thisiter != TempData.end())
		{
			return true;
		}
		return false;
	};

	//void Element::CopytoArray(double *array) const
	//{
	//	integer N = length, inc = 1;
	//	dcopy_(&N, Space, &inc, array, &inc);
	//};

	void Element::ObtainTempNames(std::string *names) const
	{
		MAP::const_iterator thisiter;
		integer idx = 0;
		for (thisiter = TempData.begin(); thisiter != TempData.end(); thisiter++, idx++)
		{
			names[idx].assign(thisiter->first);
		}
	};
}; /*end of ROPTLIB namespace*/
