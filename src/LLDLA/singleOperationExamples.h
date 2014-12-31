#include "DLAReg.h"
#include "LLDLA.h"
#include "LLDLATranspose.h"
#include "madd.h"
#include "mvmul.h"
#include "smmul.h"
#include "svmul.h"
#include "vmmul.h"
#include "vadd.h"
#include "vvdot.h"

#if DOLLDLA

RealPSet* GemmExample(Type dataType, Trans transA, Trans transB, int m, int n, int p);
RealPSet* DotExample(Type dataType, int m);
RealPSet* MVMulExample(Type dataType, bool transpose, int m, int n);
RealPSet* MAddExample(Type dataType, int m, int n);
RealPSet* SVMulRowExample(Type dataType, int m);
RealPSet* SVMulColExample(Type dataType, int m);
RealPSet* SMMulExample(Type dataType, int m, int n);
RealPSet* VMMulExample(Type dataType, int m, int n);
RealPSet* VAddExample(Type dataType, int m);

#endif // DOLLDLA
