#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int ssyrfs_(char *uplo, integer *n, integer *nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif