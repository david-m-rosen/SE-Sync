#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sgtsvx_(char *fact, char *trans, integer *n, integer *nrhs, real *dl, real *d__, real *du, real *dlf, real *df, real *duf, real *du2, integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif