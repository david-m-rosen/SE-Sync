#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

doublereal dlamch_(char *cmach);
int dlamc1_(integer *beta, integer *t, logical *rnd, logical *ieee1);
int dlamc2_(integer *beta, integer *t, logical *rnd, doublereal *eps, integer *emin, doublereal *rmin, integer *emax, doublereal *rmax);
doublereal dlamc3_(doublereal *a, doublereal *b);
int dlamc4_(integer *emin, doublereal *start, integer *base);
int dlamc5_(integer *beta, integer *p, integer *emin, logical *ieee, integer *emax, doublereal *rmax);

#ifdef __cplusplus
}
#endif