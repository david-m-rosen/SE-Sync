#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgemv_(char *trans, integer *m, integer *n, complex *alpha, complex *a, integer *lda, complex *x, integer *incx, complex *beta, complex *y, integer *incy);

#ifdef __cplusplus
}
#endif