#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int ssysvx_(char *fact, char *uplo, integer *n, integer *nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, integer *lwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif