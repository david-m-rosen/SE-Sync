#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zhpgvx_(integer *itype, char *jobz, char *range, char *uplo, integer *n, doublecomplex *ap, doublecomplex *bp, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__, integer *ldz, doublecomplex *work, doublereal *rwork, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif