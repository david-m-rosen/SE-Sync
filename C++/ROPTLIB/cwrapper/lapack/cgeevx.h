#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, complex *a, integer *lda, complex *w, complex *vl, integer *ldvl, complex *vr, integer *ldvr, integer *ilo, integer *ihi, real *scale, real *abnrm, real *rconde, real *rcondv, complex *work, integer *lwork, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif