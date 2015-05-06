#include "LLDLA.h"

#ifndef BLAS_EXAMPLES_H_
#define BLAS_EXAMPLES_H_

#if DOLLDLA

RealPSet* TRSMLTest(Type dataType, int m, int n);
RealPSet* GemmTest(Type dataType, Trans transA, Trans transB, int m, int n, int p);
RealPSet* Axpy(Type dataType, VecType vType, int m);
RealPSet* Gemv(Type dataType, bool transpose, int m, int n);

#endif //DOLLDLA

#endif // BLAS_EXAMPLES_H_
