#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zherfs_(char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#ifdef __cplusplus
}
#endif