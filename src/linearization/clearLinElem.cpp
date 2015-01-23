#include "clearLinElem.h"

void ClearLinElem::Print(IndStream &out)
{
#if DOTENSORS
  out.Indent();
  *out << m_name << ".EmptyData();\n";
#else
  throw;
#endif  
}
