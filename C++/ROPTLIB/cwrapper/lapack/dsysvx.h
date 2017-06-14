#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dsysvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *lwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif