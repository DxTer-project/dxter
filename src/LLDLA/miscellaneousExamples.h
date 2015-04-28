#include "DLAReg.h"
#include "LLDLA.h"

#if DOLLDLA

RealPSet* BasicMultiAssign(Type dataType, int m, int n);
RealPSet* LGenCompareL2(Type dataType, int m, int n);
RealPSet* LGenCompareL1(Type dataType, int m);
RealPSet* PackTest(Type dataType, int m);
RealPSet* SetToZeroTest(Type dataType, int m, int n);
RealPSet* GenSizeColSVMul(Type dataType, int m);
RealPSet* TransMVMul(Type dataType, int m, int n);
RealPSet* MVMul2(Type dataType, int m, int n, int p);
RealPSet* MAdd2(Type dataType, int m, int n);
RealPSet* VAdd2(Type dataType, int m);
RealPSet* DoubleGemm(Type dataType, Trans transA, Trans transB, int m, int n, int p, int k);

#endif // DOLLDLA
