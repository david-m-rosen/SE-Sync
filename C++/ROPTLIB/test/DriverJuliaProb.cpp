
#include "DriverJuliaProb.h"


#ifdef DRIVERJULIAPROB

using namespace ROPTLIB;

double *DriverJuliaProb(struct FunHandles *Handles, struct SolverParams *Sparams,
                     struct ManiParams *Mparams, integer inHasHHR, double *X0, integer length_X0, double *soln)
{
    // Initialization of Julia is not necessary.
    // If this code is run in C++ environment,  then calling julia requires initialization of Julia.
    // However, this C++ code is run in Julia. This implies Julia has been run when this code is called.
    // Therefore, it is not necessary to initialize Julia.
    //	jl_init(JULIA_DIR)

    // Obtain manifold and iterate structure
    Manifold *domain, **manifolds;
    Variable *initialX;
    Element **elements;
    integer *powsinterval, numoftype, numoftotal;

    if (!ParseManiParams(Mparams, manifolds, elements, powsinterval, numoftype, numoftotal))
    {
        std::cout << "Parsing ManiParams fails." << std::endl;
        exit(EXIT_FAILURE);
    }

    domain = new ProductManifold(manifolds, numoftype, powsinterval, numoftotal);

    bool HasHHR = (inHasHHR != 0);
    domain->SetHasHHR(HasHHR);

    if(Mparams->IsCheckParams != 0)
        domain->CheckParams();

    initialX = new ProductElement(elements, numoftotal, powsinterval, numoftype);
    //initialX->Print("initialX");

    // initialize the initial iterate
    if (length_X0 != 0)
    {
        double *initXptr = initialX->ObtainWriteEntireData();
        if(length_X0 != initialX->Getlength())
        {
            std::cout << "Error: The initial iterate does not have correct size: " << length_X0 << "!=" << initialX->Getlength() << "!" << std::endl;
            exit(EXIT_FAILURE);
        }
        dcopy_(&length_X0, X0, &GLOBAL::IONE, initXptr, &GLOBAL::IONE);
    }
    else
    {
        initialX->RandInManifold();
	}
	Variable *SolnX = nullptr;
	if (soln != nullptr)
	{
		SolnX = new ProductElement(elements, numoftotal, powsinterval, numoftype);
		double *SolnXptr = SolnX->ObtainWriteEntireData();
		dcopy_(&length_X0, soln, &GLOBAL::IONE, SolnXptr, &GLOBAL::IONE);
	}

    // Define the problem
    jl_function_t *func = jl_get_function(jl_main_module, Handles->fname);
    jl_function_t *gfunc = jl_get_function(jl_main_module, Handles->gfname);
    jl_function_t *hfunc = jl_get_function(jl_main_module, Handles->hfname);

    Problem *Prob = new juliaProblem(func, gfunc, hfunc);
    Prob->SetDomain(domain);

    double *result = ParseSolverParamsAndOptimizing(Sparams, Handles, Prob, initialX, SolnX);

    delete Prob;
    delete domain;
	delete initialX;
	delete SolnX;

    for (integer i = 0; i < numoftype; i++)
    {
        delete manifolds[i];
        delete elements[i];
    }
    delete[] manifolds;
    delete[] elements;
    delete[] powsinterval;

    return result;

//    It is not necessary to close Julia. The reason is the same as that in "jl_init"
//    jl_atexit_hook(0);
//    return;
};

namespace RJULIA{
    jl_function_t *isstopped = nullptr;
    /*This function defines the stopping criterion that may be used in the C++ solver*/
    bool juliaInnerStop(Variable *x, Vector *gf, double f, double ngf, double ngf0, const Problem *prob, const Solvers *solver)
    {
        jl_value_t *args[5] = {nullptr};
        jl_value_t* array_type = jl_apply_array_type(jl_float64_type, 1);
        double *xptr = x->ObtainWritePartialData();
        jl_array_t *arrx = jl_ptr_to_array_1d(array_type, xptr, x->Getlength(), 0);
        double *gfptr = gf->ObtainWritePartialData();
        jl_array_t *arrgf = jl_ptr_to_array_1d(array_type, gfptr, gf->Getlength(), 0);
        args[0] = (jl_value_t *) arrx;
        args[1] = (jl_value_t *) arrgf;
        args[2] = jl_box_float64(f);
        args[3] = jl_box_float64(ngf);
        args[4] = jl_box_float64(ngf0);

        jl_value_t *retx = jl_call(isstopped, args, 5);

        if(jl_is_int64(retx))
        {
            integer result = jl_unbox_int64(retx);
            return (result != 0);
        }

        if(jl_is_bool(retx))
        {
            return jl_unbox_bool(retx);
        }

        std::cout << "Error: Function isstopped must return an integer or a bool!" << std::endl;
        exit(EXIT_FAILURE);
    };

    jl_function_t *LinesearchInput = nullptr;
    /*This function defines the line search algorithm that may be used in the C++ solver*/
    double juliaLinesearchInput(integer iter, Variable *x1, Vector *eta1, double initialstepsize, double initialslope, const Problem *prob, const Solvers *solver)
    {
        jl_value_t *args[4] = {nullptr};
        jl_value_t* array_type = jl_apply_array_type(jl_float64_type, 1);
        double *xptr = x1->ObtainWritePartialData();
        jl_array_t *arrx = jl_ptr_to_array_1d(array_type, xptr, x1->Getlength(), 0);
        double *etaptr = eta1->ObtainWritePartialData();
        jl_array_t *arreta = jl_ptr_to_array_1d(array_type, etaptr, eta1->Getlength(), 0);
        args[0] = (jl_value_t *) arrx;
        args[1] = (jl_value_t *) arreta;
        args[2] = jl_box_float64(initialstepsize);
		args[3] = jl_box_float64(initialslope);
		args[4] = jl_box_int64(iter);

        jl_value_t *retx = jl_call(LinesearchInput, args, 5);

        if(jl_is_float64(retx))
        {
            double result = jl_unbox_float64(retx);
            return result;
        }

        std::cout << "Error: Function isstopped must return a number of double precision!" << std::endl;
        exit(EXIT_FAILURE);
    }
};

double *ParseSolverParamsAndOptimizing(struct SolverParams *Sparams, struct FunHandles *Handles, Problem *Prob, Variable *initialX, Variable *SolnX)
{
    PARAMSMAP params;

//    Add parameters
    if(Sparams->Stop_Criterion != -1)
        params.insert(std::pair<std::string, double>("Stop_Criterion", Sparams->Stop_Criterion));
    if(Sparams->Tolerance != -1)
        params.insert(std::pair<std::string, double>("Tolerance", Sparams->Tolerance));
    if(Sparams->Diffx != -1)
        params.insert(std::pair<std::string, double>("Diffx", Sparams->Diffx));
    if(Sparams->NumExtraGF != -1)
        params.insert(std::pair<std::string, double>("NumExtraGF", Sparams->NumExtraGF));
    if(Sparams->TimeBound != -1)
        params.insert(std::pair<std::string, double>("TimeBound", Sparams->TimeBound));
    if(Sparams->Min_Iteration != -1)
        params.insert(std::pair<std::string, double>("Min_Iteration", Sparams->Min_Iteration));
    if(Sparams->Max_Iteration != -1)
        params.insert(std::pair<std::string, double>("Max_Iteration", Sparams->Max_Iteration));
    if(Sparams->OutputGap != -1)
        params.insert(std::pair<std::string, double>("OutputGap", Sparams->OutputGap));
    if(Sparams->DEBUG != -1)
        params.insert(std::pair<std::string, double>("DEBUG", Sparams->DEBUG));
    if(Sparams->isconvex != -1)
        params.insert(std::pair<std::string, double>("isconvex", Sparams->isconvex));
    if(Sparams->nu != -1)
        params.insert(std::pair<std::string, double>("nu", Sparams->nu));
    if(Sparams->mu != -1)
        params.insert(std::pair<std::string, double>("mu", Sparams->mu));
    if(Sparams->LengthSY != -1)
        params.insert(std::pair<std::string, double>("LengthSY", Sparams->LengthSY));
    if(Sparams->lambdaLower != -1)
        params.insert(std::pair<std::string, double>("lambdaLower", Sparams->lambdaLower));
    if(Sparams->lambdaUpper != -1)
        params.insert(std::pair<std::string, double>("lambdaUpper", Sparams->lambdaUpper));
    if(Sparams->LineSearch_LS != -1)
        params.insert(std::pair<std::string, double>("LineSearch_LS", Sparams->LineSearch_LS));
    if(Sparams->IsPureLSInput != -1)
        params.insert(std::pair<std::string, double>("IsPureLSInput", Sparams->IsPureLSInput));
    if(Sparams->LS_alpha != -1)
        params.insert(std::pair<std::string, double>("LS_alpha", Sparams->LS_alpha));
    if(Sparams->LS_beta != -1)
        params.insert(std::pair<std::string, double>("LS_beta", Sparams->LS_beta));
    if(Sparams->Minstepsize != -1)
        params.insert(std::pair<std::string, double>("Minstepsize", Sparams->Minstepsize));
    if(Sparams->Maxstepsize != -1)
        params.insert(std::pair<std::string, double>("Maxstepsize", Sparams->Maxstepsize));
    if(Sparams->LS_ratio1 != -1)
        params.insert(std::pair<std::string, double>("LS_ratio1", Sparams->LS_ratio1));
    if(Sparams->LS_ratio2 != -1)
        params.insert(std::pair<std::string, double>("LS_ratio2", Sparams->LS_ratio2));
    if(Sparams->Initstepsize != -1)
        params.insert(std::pair<std::string, double>("Initstepsize", Sparams->Initstepsize));
    if(Sparams->Accuracy != -1)
        params.insert(std::pair<std::string, double>("Accuracy", Sparams->Accuracy));
    if(Sparams->Finalstepsize != -1)
        params.insert(std::pair<std::string, double>("Finalstepsize", Sparams->Finalstepsize));
    if(Sparams->Num_pre_funs != -1)
        params.insert(std::pair<std::string, double>("Num_pre_funs", Sparams->Num_pre_funs));
    if(Sparams->InitSteptype != -1)
        params.insert(std::pair<std::string, double>("InitSteptype", Sparams->InitSteptype));
    if(Sparams->Acceptence_Rho != -1)
        params.insert(std::pair<std::string, double>("Acceptence_Rho", Sparams->Acceptence_Rho));
    if(Sparams->Shrinked_tau != -1)
        params.insert(std::pair<std::string, double>("Shrinked_tau", Sparams->Shrinked_tau));
    if(Sparams->Magnified_tau != -1)
        params.insert(std::pair<std::string, double>("Magnified_tau", Sparams->Magnified_tau));
    if(Sparams->minimum_Delta != -1)
        params.insert(std::pair<std::string, double>("minimum_Delta", Sparams->minimum_Delta));
    if(Sparams->maximum_Delta != -1)
        params.insert(std::pair<std::string, double>("maximum_Delta", Sparams->maximum_Delta));
    if(Sparams->useRand != -1)
        params.insert(std::pair<std::string, double>("useRand", Sparams->useRand));
    if(Sparams->Max_Inner_Iter != -1)
        params.insert(std::pair<std::string, double>("Max_Inner_Iter", Sparams->Max_Inner_Iter));
    if(Sparams->Min_Inner_Iter != -1)
        params.insert(std::pair<std::string, double>("Min_Inner_Iter", Sparams->Min_Inner_Iter));
    if(Sparams->theta != -1)
        params.insert(std::pair<std::string, double>("theta", Sparams->theta));
    if(Sparams->kappa != -1)
        params.insert(std::pair<std::string, double>("kappa", Sparams->kappa));
    if(Sparams->initial_Delta != -1)
        params.insert(std::pair<std::string, double>("initial_Delta", Sparams->initial_Delta));
    if(Sparams->Eps != -1)
        params.insert(std::pair<std::string, double>("Eps", Sparams->Eps));
    if(Sparams->Theta_eps != -1)
        params.insert(std::pair<std::string, double>("Theta_eps", Sparams->Theta_eps));
    if(Sparams->Min_Eps != -1)
        params.insert(std::pair<std::string, double>("Min_Eps", Sparams->Min_Eps));
    if(Sparams->Del != -1)
        params.insert(std::pair<std::string, double>("Del", Sparams->Del));
    if(Sparams->Theta_del != -1)
        params.insert(std::pair<std::string, double>("Theta_del", Sparams->Theta_del));

    std::string stdmethodname = Sparams->name;
    Solvers *solver;
    if (stdmethodname == "RSD")
    {
        solver = new RSD(Prob, initialX, SolnX);
    }
    else
    if (stdmethodname == "RNewton")
    {
		solver = new RNewton(Prob, initialX, SolnX);
    }
    else
    if (stdmethodname == "RCG")
    {
		solver = new RCG(Prob, initialX, SolnX);
    }
    else
    if (stdmethodname == "RBroydenFamily")
    {
		solver = new RBroydenFamily(Prob, initialX, nullptr, SolnX);
    }
    else
    if (stdmethodname == "RWRBFGS")
    {
		solver = new RWRBFGS(Prob, initialX, nullptr, SolnX);
    }
    else
    if (stdmethodname == "RBFGS")
    {
		solver = new RBFGS(Prob, initialX, nullptr, SolnX);
    }
    else
    if (stdmethodname == "RBFGSLPSub")
    {
		solver = new RBFGSLPSub(Prob, initialX, nullptr, SolnX);
    }
    else
    if (stdmethodname == "LRBFGSLPSub")
    {
		solver = new LRBFGSLPSub(Prob, initialX, SolnX);
    }
    else
    if (stdmethodname == "RGS")
    {
		solver = new RGS(Prob, initialX, SolnX);
    }
    else
    if (stdmethodname == "LRBFGS")
    {
		solver = new LRBFGS(Prob, initialX, SolnX);
    }
    else
    if (stdmethodname == "RTRSD")
    {
		solver = new RTRSD(Prob, initialX, SolnX);
    }
    else
    if (stdmethodname == "RTRNewton")
    {
		solver = new RTRNewton(Prob, initialX, SolnX);
    }
    else
    if (stdmethodname == "RTRSR1")
    {
		solver = new RTRSR1(Prob, initialX, nullptr, SolnX);
    }
    else
    if (stdmethodname == "LRTRSR1")
    {
		solver = new LRTRSR1(Prob, initialX, SolnX);
    }
    else
    {
        std::cout << "Warning: Unrecognized solver: " << stdmethodname << ". Use LRBFGS instead!" << std::endl;
        solver = new LRBFGS(Prob, initialX);
    }
    solver->SetParams(params);

    if(strlen(Handles->isstopped) > 0)
    {
        RJULIA::isstopped = jl_get_function(jl_main_module, Handles->isstopped);
        solver->StopPtr = &RJULIA::juliaInnerStop;
    }
    if(strlen(Handles->LinesearchInput) > 0)
    {
        RJULIA::LinesearchInput = jl_get_function(jl_main_module, Handles->LinesearchInput);
        SolversLS *solverLS = dynamic_cast<SolversLS *> (solver);
        if (solverLS != nullptr)
            solverLS->LinesearchInput = &RJULIA::juliaLinesearchInput;
    }

    if(Sparams->IsCheckParams != 0)
        solver->CheckParams();
    solver->Run();

    integer lengthSeries = solver->GetlengthSeries();
    integer lengthresult = 1 + initialX->Getlength() + 11 + 4 * lengthSeries;
    double *result = new double[lengthresult];
    result[0] = lengthresult;
    const double *xoptptr = solver->GetXopt()->ObtainReadData();
    integer lengthx = initialX->Getlength();
    dcopy_(&lengthx, const_cast<double *>(xoptptr), &GLOBAL::IONE, result + 1, &GLOBAL::IONE);
    result[lengthx + 1] = static_cast<double> (solver->Getfinalfun());
    result[lengthx + 2] = static_cast<double> (solver->Getnormgf());
    result[lengthx + 3] = static_cast<double> (solver->Getnormgfgf0());
    result[lengthx + 4] = static_cast<double> (solver->GetIter());
    result[lengthx + 5] = static_cast<double> (solver->Getnf());
    result[lengthx + 6] = static_cast<double> (solver->Getng());
    result[lengthx + 7] = static_cast<double> (solver->GetnR());
    result[lengthx + 8] = static_cast<double> (solver->GetnV());
    result[lengthx + 9] = static_cast<double> (solver->GetnVp());
    result[lengthx + 10] = static_cast<double> (solver->GetnH());
    result[lengthx + 11] = static_cast<double> (solver->GetComTime());
    for (integer i = 0; i < lengthSeries; i++)
        result[lengthx + 12 + i] = solver->GetfunSeries()[i];
    for (integer i = 0; i < lengthSeries; i++)
		result[lengthx + 12 + lengthSeries + i] = solver->GetgradSeries()[i];
	for (integer i = 0; i < lengthSeries; i++)
		result[lengthx + 12 + 2 * lengthSeries + i] = solver->GettimeSeries()[i];
	for (integer i = 0; i < lengthSeries; i++)
		result[lengthx + 12 + 3 * lengthSeries + i] = solver->GetdistSeries()[i];

    if(Sparams->IsCheckGradHess != 0)
    {
        Prob->CheckGradHessian(initialX);
        // Check gradient and Hessian
        const Variable *xopt = solver->GetXopt();
        Variable *xoptcopy = xopt->ConstructEmpty();
        xopt->CopyTo(xoptcopy);
        xoptcopy->RemoveAllFromTempData();
        Prob->CheckGradHessian(xoptcopy);
        delete xoptcopy;
    }
    delete solver;
    return result;
};

bool ParseManiParams(struct ManiParams *Mparams, Manifold **&manifolds, Element **&elements,
                     integer *&powsinterval, integer &numoftype, integer &numoftotal)
{
    // Parse ManiParams
    numoftype = Mparams->numoftypes;

    powsinterval = new integer[numoftype + 1];
    const char *name = nullptr;
    manifolds = new Manifold *[numoftype];
    integer n, p, m, Params;
    powsinterval[0] = 0;

    for (integer i = 0; i < numoftype; i++)
        powsinterval[i + 1] = powsinterval[i] + Mparams->numofmani[i];

    numoftotal = powsinterval[numoftype];
    elements = new Element *[numoftotal];
    PARAMSMAP params;
    for (integer i = 0; i < numoftype; i++)
    {
        name = Mparams->name[i];
        n = Mparams->n[i];
        m = Mparams->m[i];
        p = Mparams->p[i];
        Params = Mparams->paramset[i];
        if(Params != -1)
            params[static_cast<std::string> ("ParamSet")] = Params;

        manifolds[i] = GetAManifold(name, n, m, p);
        manifolds[i]->SetParams(params);

        for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
        {
            elements[j] = GetAnElement(name, n, m, p);
        }
        if (manifolds[i] == nullptr || elements[i] == nullptr)
        {
            return false;
        }
    }
    return true;
};

Manifold *GetAManifold(const char *name, integer n, integer m, integer p)
{
    if (strcmp(name, "Euclidean") == 0)
    {
        return new Euclidean(n, m);
    }
    else
    if (strcmp(name, "Sphere") == 0)
    {
        return new Sphere(n);
    }
    else
    if (strcmp(name, "Stiefel") == 0)
    {
        return new Stiefel(n, p);
    }
    else
    if (strcmp(name, "Oblique") == 0)
    {
        return new Oblique(n, m);
    }
    else
    if (strcmp(name, "LowRank") == 0)
    {
        return new LowRank(n, m, p);
    }
    else
    if (strcmp(name, "OrthGroup") == 0)
    {
        return new OrthGroup(n);
    }
    else
    if (strcmp(name, "L2Sphere") == 0)
    {
        return new L2Sphere(n);
    }
    else
    if (strcmp(name, "SPDManifold") == 0)
    {
        return new SPDManifold(n);
    }
    else
    if (strcmp(name, "CpxNStQOrth") == 0)
    {
        return new CpxNStQOrth(n, p);
    }
    else
    if (strcmp(name, "Grassmann") == 0)
    {
        return new Grassmann(n, p);
    }
    else
    if (strcmp(name, "EucPositive") == 0)
    {
        return new EucPositive(n, m);
    }
    else
    if (strcmp(name, "SPDTensor") == 0)
    {
        return new SPDTensor(n, m);
    }
    else
    {
        std::cout << "Manifold: " << name << " does not implemented in this library!" << std::endl;
        return nullptr;
    }
};

Element *GetAnElement(const char *name, integer n, integer m, integer p)
{
    if (strcmp(name, "Euclidean") == 0)
    {
        return new EucVariable(n, m);
    }
    else
    if (strcmp(name, "Sphere") == 0)
    {
        return new SphereVariable(n);
    }
    else
    if (strcmp(name, "Stiefel") == 0)
    {
        return new StieVariable(n, p);
    }
    else
    if (strcmp(name, "Oblique") == 0)
    {
        return new ObliqueVariable(n, m);
    }
    else
    if (strcmp(name, "LowRank") == 0)
    {
        return new LowRankVariable(n, m, p);
    }
    else
    if (strcmp(name, "OrthGroup") == 0)
    {
        return new OrthGroupVariable(n);
    }
    else
    if (strcmp(name, "L2Sphere") == 0)
    {
        return new L2SphereVariable(n);
    }
    else
    if (strcmp(name, "SPDManifold") == 0)
    {
        return new SPDVariable(n);
    }
    else
    if (strcmp(name, "CpxNStQOrth") == 0)
    {
        return new CSOVariable(n, p);
    }
    else
    if (strcmp(name, "Grassmann") == 0)
    {
        return new GrassVariable(n, p);
    }
    else
    if (strcmp(name, "EucPositive") == 0)
    {
        return new EucPosVariable(n, m);
    }
    else
    if (strcmp(name, "SPDTensor") == 0)
    {
        return new SPDTVariable(n, m);
    }
    else
    {
        std::cout << "Element: " << name << " does not implemented in this library!" << std::endl;
        return nullptr;
    }
};

#endif // end of DRIVERJULIAPROB
