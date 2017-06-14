#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cspsvx_(char *fact, char *uplo, integer *n, integer *nrhs, complex *ap, complex *afp, integer *ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif