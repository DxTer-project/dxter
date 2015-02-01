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

RealPSet* GemmTest(Type dataType, Trans transA, Trans transB, int m, int n, int p);
RealPSet* DotTest(Type dataType, int m);
RealPSet* MVMulTest(Type dataType, bool transpose, int m, int n);
RealPSet* MAddTest(Type dataType, int m, int n);
RealPSet* SVMulTest(Type dataType, VecType vecType, int m);
RealPSet* SMMulTest(Type dataType, int m, int n);
RealPSet* VMMulTest(Type dataType, int m, int n);
RealPSet* VAddTest(Type dataType, VecType vecType, int m);

#endif // DOLLDLA
