#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int ztprfs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#ifdef __cplusplus
}
#endif