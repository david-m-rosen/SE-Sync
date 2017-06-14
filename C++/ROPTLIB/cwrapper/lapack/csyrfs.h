#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int csyrfs_(char *uplo, integer *n, integer *nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif