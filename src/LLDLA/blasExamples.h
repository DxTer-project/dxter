#include "LLDLA.h"

#ifndef BLAS_EXAMPLES_H_
#define BLAS_EXAMPLES_H_

#if DOLLDLA

RealPSet* Axpy(Type dataType, VecType vType, int m);
RealPSet* Gemv(Type dataType, bool transpose, int m, int n);

#endif //DOLLDLA

#endif // BLAS_EXAMPLES_H_
