#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgbmv_(char *trans, integer *m, integer *n, integer *kl, integer *ku, complex *alpha, complex *a, integer *lda, complex *x, integer *incx, complex *beta, complex *y, integer *incy);

#ifdef __cplusplus
}
#endif