#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dgtsvx_(char *fact, char *trans, integer *n, integer *nrhs, doublereal *dl, doublereal *d__, doublereal *du, doublereal *dlf, doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif