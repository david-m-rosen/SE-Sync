#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dtprfs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif