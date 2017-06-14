#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zgemv_(char *trans, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy);

#ifdef __cplusplus
}
#endif