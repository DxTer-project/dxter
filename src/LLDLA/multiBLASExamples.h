#include "DLAReg.h"
#include "LLDLA.h"
#include "LLDLATranspose.h"
#include "madd.h"
#include "mvmul.h"
#include "smmul.h"
#include "svmul.h"
#include "vadd.h"
#include "vvdot.h"

#if DOLLDLA

RealPSet* Gesummv(Type dataType, int m, int n);
RealPSet* Gemam(Type dataType, int m, int n, int p);
RealPSet* Gemv2(Type dataType, int m, int n, int k);
RealPSet* VMVMulExample(Type dataType, int m, int n);

#endif //DOLLDLA
