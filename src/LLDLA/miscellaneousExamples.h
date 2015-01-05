#include "DLAReg.h"
#include "LLDLA.h"
#include "madd.h"
#include "mvmul.h"
#include "partition.h"
#include "recombine.h"
#include "svmul.h"
#include "vadd.h"
#include "vvdot.h"

#if DOLLDLA

RealPSet* GenSizeColSVMul(Type dataType, int m);
RealPSet* TransMVMulExample(Type dataType, int m, int n);
RealPSet* MVMul2Example(Type dataType, int m, int n, int p);
RealPSet* MAdd2Example(Type dataType, int m, int n);
RealPSet* VAdd2Example(Type dataType, int m);
RealPSet* DoubleGemmExample(Type dataType, Trans transA, Trans transB, int m, int n, int p, int k);

#endif // DOLLDLA
