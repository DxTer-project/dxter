#include "DLAReg.h"
#include "LLDLA.h"
#include "LLDLATranspose.h"
#include "mvmul.h"
#include "svmul.h"
#include "vadd.h"

#if DOLLDLA

RealPSet* GemvExample(Type dataType, bool transpose, int m, int n);

#endif //DOLLDLA
