#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dgbsvx_(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, integer *ipiv, char *equed, doublereal *r__, doublereal *c__, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif