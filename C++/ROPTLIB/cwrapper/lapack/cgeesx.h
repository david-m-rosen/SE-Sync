#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgeesx_(char *jobvs, char *sort, L_fp select, char *sense, integer *n, complex *a, integer *lda, integer *sdim, complex *w, complex *vs, integer *ldvs, real *rconde, real *rcondv, complex *work, integer *lwork, real *rwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif