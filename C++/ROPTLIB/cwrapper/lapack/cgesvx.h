#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgesvx_(char *fact, char *trans, integer *n, integer *nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *ipiv, char *equed, real *r__, real *c__, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif