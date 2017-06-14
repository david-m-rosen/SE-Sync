
#include "juliaProblem.h"

#ifdef DRIVERJULIAPROB

/*Define the namespace*/
namespace ROPTLIB{

    juliaProblem::juliaProblem(jl_function_t *inf, jl_function_t *ingf, jl_function_t *inHess)
    {
        jl_f = inf;
        jl_gf = ingf;
        jl_Hess = inHess;
	};

    juliaProblem::~juliaProblem()
	{
	};

    double juliaProblem::f(Variable *x) const
    {
//        x->Print("cpp f x");//---
        jl_value_t* array_type = jl_apply_array_type(jl_float64_type, 1);
        double *xptr = x->ObtainWritePartialData();
        jl_array_t *arrx = jl_ptr_to_array_1d(array_type, xptr, x->Getlength(), 0);

        jl_array_t *arrtmp = nullptr;
        if(x->TempDataExist(("Tmp")))
        {
            const SharedSpace *Tmp = x->ObtainReadTempData("Tmp");
            const double *tmpptr = Tmp->ObtainReadData();
            arrtmp = jl_ptr_to_array_1d(array_type, const_cast<double *> (tmpptr), Tmp->Getlength(), 0);
        } else
        {
            arrtmp = jl_ptr_to_array_1d(array_type, nullptr, 0, 0);
        }

        jl_value_t *retresult = jl_call2(jl_f, (jl_value_t *) arrx, (jl_value_t *) arrtmp);
        jl_get_nth_field(retresult, 0);
        jl_value_t *fx = jl_get_nth_field(retresult, 0);
        jl_array_t *outtmp = (jl_array_t *) jl_get_nth_field(retresult, 1);

        integer outtmplen = jl_array_len(outtmp);
        SharedSpace *sharedouttmp = new SharedSpace(1, outtmplen);
        double *outtmpptr = sharedouttmp->ObtainWriteEntireData();
        dcopy_(&outtmplen, (double*)jl_array_data(outtmp), &GLOBAL::IONE, outtmpptr, &GLOBAL::IONE);
//        sharedouttmp->Print("cpp f tmp:");//----
        x->RemoveFromTempData("Tmp");
        x->AddToTempData("Tmp", sharedouttmp);

        if(jl_is_float64(fx))
        {
            double result = jl_unbox_float64(fx);
//            std::cout << "cpp f fx:" << result << std::endl;//-----
            return result;
        }
        std::cout << "Error: The objectve function must return a number of double precision!" << std::endl;
        exit(EXIT_FAILURE);
	};

    void juliaProblem::EucGrad(Variable *x, Vector *egf) const
    {
//        x->Print("cpp gf x");//---
        jl_value_t* array_type = jl_apply_array_type(jl_float64_type, 1);
        double *xptr = x->ObtainWritePartialData();
        jl_array_t *arrx = jl_ptr_to_array_1d(array_type, xptr, x->Getlength(), 0);

        jl_array_t *arrtmp = nullptr;
        if(x->TempDataExist(("Tmp")))
        {
            const SharedSpace *Tmp = x->ObtainReadTempData("Tmp");
//            Tmp->Print("cpp gf inTmp");//---
            const double *tmpptr = Tmp->ObtainReadData();
            arrtmp = jl_ptr_to_array_1d(array_type, const_cast<double *> (tmpptr), Tmp->Getlength(), 0);
        } else
        {
            arrtmp = jl_ptr_to_array_1d(array_type, nullptr, 0, 0);
        }

        jl_value_t *retresult = jl_call2(jl_gf, (jl_value_t *) arrx, (jl_value_t *) arrtmp);
        jl_array_t *jl_egf = (jl_array_t *) jl_get_nth_field(retresult, 0);
        jl_array_t *outtmp = (jl_array_t *) jl_get_nth_field(retresult, 1);

        if(jl_array_len(jl_egf) != egf->Getlength())
        {
            std::cout << "error: the size of the Euclidean gradient is not correct!" << std::endl;
            exit(EXIT_FAILURE);
        }

        integer egflen = egf->Getlength();
        double *egfptr = egf->ObtainWriteEntireData();
        dcopy_(&egflen, (double*)jl_array_data(jl_egf), &GLOBAL::IONE, egfptr, &GLOBAL::IONE);
//        egf->Print("cpp gf egf:");//--

        integer outtmplen = jl_array_len(outtmp);
        if(outtmplen != 0)
        {
            SharedSpace *sharedouttmp = new SharedSpace(1, outtmplen);
            double *outtmpptr = sharedouttmp->ObtainWriteEntireData();
            dcopy_(&outtmplen, (double*)jl_array_data(outtmp), &GLOBAL::IONE, outtmpptr, &GLOBAL::IONE);
            x->RemoveFromTempData("Tmp");
            x->AddToTempData("Tmp", sharedouttmp);
        }
	};

    void juliaProblem::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
    {
//        x->Print("cpp hf x");//---
//        etax->Print("cpp hf etax");//---
        jl_value_t* array_type = jl_apply_array_type(jl_float64_type, 1);
        double *xptr = x->ObtainWritePartialData();
        jl_array_t *arrx = jl_ptr_to_array_1d(array_type, xptr, x->Getlength(), 0);
        double *etaxptr = etax->ObtainWritePartialData();
        jl_array_t *arretax = jl_ptr_to_array_1d(array_type, etaxptr, etax->Getlength(), 0);

        jl_array_t *arrtmp = nullptr;
        if(x->TempDataExist(("Tmp")))
        {
            const SharedSpace *Tmp = x->ObtainReadTempData("Tmp");
//            Tmp->Print("cpp hf inTmp");//---
            const double *tmpptr = Tmp->ObtainReadData();
            arrtmp = jl_ptr_to_array_1d(array_type, const_cast<double *> (tmpptr), Tmp->Getlength(), 0);
        } else
        {
            arrtmp = jl_ptr_to_array_1d(array_type, nullptr, 0, 0);
        }

        jl_value_t *retresult = jl_call3(jl_Hess, (jl_value_t *) arrx, (jl_value_t *) arrtmp, (jl_value_t *) arretax);
        jl_array_t *jl_exix = (jl_array_t *) jl_get_nth_field(retresult, 0);
        jl_array_t *outtmp = (jl_array_t *) jl_get_nth_field(retresult, 1);

        if(jl_array_len(jl_exix) != etax->Getlength())
        {
            std::cout << "error: the size of the action of the Hessian is not correct!" << std::endl;
            exit(EXIT_FAILURE);
        }

        integer exixlen = exix->Getlength();
        double *exixptr = exix->ObtainWriteEntireData();
        dcopy_(&exixlen, (double*)jl_array_data(jl_exix), &GLOBAL::IONE, exixptr, &GLOBAL::IONE);
//        exix->Print("cpp hf exix:");//---

        integer outtmplen = jl_array_len(outtmp);
        if(outtmplen != 0)
        {
            SharedSpace *sharedouttmp = new SharedSpace(1, outtmplen);
            double *outtmpptr = sharedouttmp->ObtainWriteEntireData();
            dcopy_(&outtmplen, (double*)jl_array_data(outtmp), &GLOBAL::IONE, outtmpptr, &GLOBAL::IONE);
            x->RemoveFromTempData("Tmp");
            x->AddToTempData("Tmp", sharedouttmp);
        }
	};

}; /*end of ROPTLIB namespace*/

#endif // end of DRIVERJULIAPROB
