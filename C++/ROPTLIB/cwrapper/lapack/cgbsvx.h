#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgbsvx_(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, complex *ab, integer *ldab, complex *afb, integer *ldafb, integer *ipiv, char *equed, real *r__, real *c__, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif