#include "DLAReg.h"
#include "LLDLA.h"

#if DOLLDLA

RealPSet* Gesummv(Type dataType, int m, int n);
RealPSet* Gemam(Type dataType, int m, int n, int p);
RealPSet* Gemv2(Type dataType, int m, int n, int k);
RealPSet* VMVMul(Type dataType, int m, int n);

#endif //DOLLDLA
