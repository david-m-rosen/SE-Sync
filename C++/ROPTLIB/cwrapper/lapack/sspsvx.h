#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sspsvx_(char *fact, char *uplo, integer *n, integer *nrhs, real *ap, real *afp, integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif