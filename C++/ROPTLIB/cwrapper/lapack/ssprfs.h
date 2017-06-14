#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int ssprfs_(char *uplo, integer *n, integer *nrhs, real *ap, real *afp, integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif