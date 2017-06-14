#include "Manifolds/ProductElement.h"

/*Define the namespace*/
namespace ROPTLIB{

	ProductElement::ProductElement(Element **inelements, integer innumofelements, integer *inpowsinterval, integer innumoftypes)
	{
		ProductElementInitialization(inelements, innumofelements, inpowsinterval, innumoftypes);
	};

	void ProductElement::ProductElementInitialization(Element **inelements, integer innumofelements, integer *inpowsinterval, integer innumoftypes)
	{
		numoftypes = innumoftypes;
		powsinterval = new integer[numoftypes + 1];
		for (integer i = 0; i < numoftypes + 1; i++)
		{
			powsinterval[i] = inpowsinterval[i];
		}
		numofelements = innumofelements;
		elements = new Element *[numofelements];

		length = 0;
		for (integer i = 0; i < numofelements; i++)
		{
			length += inelements[i]->Getlength();
		}

		size = new integer; *size = length;
		ls = 1;
		sharedtimes = nullptr;
		Space = nullptr;

		for (integer i = 0; i < numoftypes; i++)
		{
			integer ils = inelements[powsinterval[i]]->Getls();
			integer *isize = new integer[ils];
			for (integer j = 0; j < ils; j++)
			{
				isize[j] = inelements[powsinterval[i]]->Getsize()[j];
			}

			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				elements[j] = inelements[j]->ConstructEmpty();
				delete[] elements[j]->Getsize();
				elements[j]->SetByParams(isize, ils, inelements[j]->Getlength(), nullptr, nullptr);
			}
		}
	};

	ProductElement::ProductElement(integer numberoftypes, ...)
	{
		numoftypes = numberoftypes;
		powsinterval = new integer[numoftypes + 1];

		Element **types = new Element *[numoftypes];

		va_list argptr;
		va_start(argptr, numberoftypes);
		powsinterval[0] = 0;
		length = 0;
		for (integer i = 0; i < numoftypes; i++)
		{
			types[i] = va_arg(argptr, Element *);
			powsinterval[i + 1] = powsinterval[i] + va_arg(argptr, integer);
			length += types[i]->Getlength() * (powsinterval[i + 1] - powsinterval[i]);
		}
		va_end(argptr);

		numofelements = powsinterval[numoftypes];
		elements = new Element *[numofelements];

		ls = 1;
		size = new integer; *size = length;
		Space = nullptr;
		sharedtimes = nullptr;
		integer *isize = nullptr;
		integer *isharedtimes = nullptr;

		for (integer i = 0; i < numoftypes; i++)
		{
			isize = new integer[types[i]->Getls()];
			for (integer j = 0; j < types[i]->Getls(); j++)
			{
				isize[j] = types[i]->Getsize()[j];
			}

			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				elements[j] = types[i]->ConstructEmpty();
				delete[] elements[j]->Getsize();
				elements[j]->SetByParams(isize, types[i]->Getls(), types[i]->Getlength(), nullptr, nullptr);
			}
		}

		delete[] types;
	};

	ProductElement::ProductElement(void)
	{
		elements = nullptr;
		numofelements = 0;

		size = nullptr;
		ls = 0;
		length = 0;
		sharedtimes = nullptr;
		Space = nullptr;
	};

	ProductElement::~ProductElement(void)
	{
		for (integer i = 0; i < numoftypes; i++)
		{
			if (elements[powsinterval[i]]->GetSharedTimes() != nullptr)
			{
				delete elements[powsinterval[i]]->GetSharedTimes();
			}
			if (elements[powsinterval[i]]->Getsize() != nullptr)
				delete[] elements[powsinterval[i]]->Getsize();
		}
		//for (integer i = 0; i < numofelements; i++)
		//{
		//	if (elements[i]->GetSharedTimes() != nullptr)
		//	{
		//		delete elements[i]->GetSharedTimes();
		//	}
		//	if (elements[i]->Getsize() != nullptr)
		//		delete[] elements[i]->Getsize();
		//	elements[i]->RemoveAllFromTempData();
		//}

		delete[] powsinterval;
		for (integer i = 0; i < numofelements; i++)
		{
			if (elements[i] != nullptr)
			{
				elements[i]->DeleteBySettingNull();
				delete elements[i];
			}
		}
		delete[] elements;
	};

	ProductElement *ProductElement::ConstructEmpty(void) const
	{
		return new ProductElement(elements, numofelements, powsinterval, numoftypes);
	};

	void ProductElement::CopyTo(Element *eta) const
	{
		if (this == eta)
			return;
		Element::CopyTo(eta);
		ProductElement *Prodeta = dynamic_cast<ProductElement *> (eta);
		Prodeta->ResetMemoryofElementsAndSpace();
	};

	void ProductElement::RandUnform(double start, double end)
	{
		ObtainWriteEntireData();
		for (integer i = 0; i < numofelements; i++)
			elements[i]->RandUnform(start, end);
	};

	void ProductElement::RandGaussian(double mean, double variance)
	{
		ObtainWriteEntireData();
		for (integer i = 0; i < numofelements; i++)
			elements[i]->RandGaussian(mean, variance);
	};

	void ProductElement::Print(const char *name, bool isonlymain) const
	{
		printf("=================Product data: %s=========================\n", name);
		for (integer i = 0; i < numofelements; i++)
		{
			//std::string s = std::to_string(i);
			std::stringstream strStream;
			strStream << "number " << i << " manifold";
			std::string str = strStream.str();
			elements[i]->Print(str.c_str(), isonlymain);
		}
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

	void ProductElement::RandInManifold(void)
	{
		this->NewMemoryOnWrite();
		for (integer i = 0; i < numoftypes; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				elements[j]->RandInManifold();
			}
		}
	};

	double *ProductElement::ObtainWriteEntireData(void)
	{
		NewMemoryOnWrite();
		return Space;
	};

	double *ProductElement::ObtainWritePartialData(void)
	{
		CopyOnWrite();
		return Space;
	};

	void ProductElement::CheckMemory(const char *info) const
	{
		for (integer i = 0; i < numoftypes; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				if (Space > elements[j]->GetSpace() || Space + length <= elements[j]->GetSpace())
				{
					printf("Error: Memory of space in elements does not use that of Space in ProductElement in %s!\n", info);
				}
			}
		}
	};

	void ProductElement::NewMemoryOnWrite(void)
	{
		Element::NewMemoryOnWrite();
		ResetMemoryofElementsAndSpace();
	};

	void ProductElement::CopyOnWrite(void)
	{
		Element::CopyOnWrite();
		ResetMemoryofElementsAndSpace();
	};

	void ProductElement::ResetMemoryofElementsAndSpace(void)
	{
		integer shift = 0;
		integer *isharedtimes = nullptr;
		ProductElement *tmp = nullptr;

		for (integer i = 0; i < numoftypes; i++)
		{
			isharedtimes = new integer;
			*isharedtimes = 1;
			if (elements[powsinterval[i]]->GetSharedTimes() != nullptr)
				delete elements[powsinterval[i]]->GetSharedTimes();
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				elements[j]->SetByParams(const_cast<integer *> (elements[j]->Getsize()), elements[j]->Getls(), elements[j]->Getlength(),
					isharedtimes, Space + shift);
				tmp = dynamic_cast<ProductElement *> (elements[j]);
				if (tmp != nullptr)
				{
					tmp->ResetMemoryofElementsAndSpace();
				}
				shift += elements[j]->Getlength();
			}
		}
	};
}; /*end of ROPTLIB namespace*/
