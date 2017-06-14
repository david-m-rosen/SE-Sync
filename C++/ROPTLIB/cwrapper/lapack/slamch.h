#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

E_f slamch_(char *cmach);
int slamc1_(integer *beta, integer *t, logical *rnd, logical *ieee1);
int slamc2_(integer *beta, integer *t, logical *rnd, real *eps, integer *emin, real *rmin, integer *emax, real *rmax);
E_f slamc3_(real *a, real *b);
int slamc4_(integer *emin, real *start, integer *base);
int slamc5_(integer *beta, integer *p, integer *emin, logical *ieee, integer *emax, real *rmax);

#ifdef __cplusplus
}
#endif