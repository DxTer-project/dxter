#include "DLAReg.h"
#include "LLDLA.h"
#include "mvmul.h"
#include "svmul.h"
#include "vadd.h"
#include "vvdot.h"

#if DOLLDLA

//RealPSet* MAddGemm(Type dataType, int m, int n, int k);
RealPSet* Gemv2(Type dataType, int m, int n, int k);
RealPSet* VMVMulExample(Type dataType, int m, int n);

#endif //DOLLDLA
