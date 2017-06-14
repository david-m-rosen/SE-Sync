#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cporfs_(char *uplo, integer *n, integer *nrhs, complex *a, integer *lda, complex *af, integer *ldaf, complex *b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif