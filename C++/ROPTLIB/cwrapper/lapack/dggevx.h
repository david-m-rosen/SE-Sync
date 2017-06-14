#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale, doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *rcondv, doublereal *work, integer *lwork, integer *iwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif