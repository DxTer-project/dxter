#include "rqoBasis.h"

#if DORQO

DataTypeInfo& DataTypeInfo::operator=(const DataTypeInfo &rhs)
{
  m_fields = rhs.m_fields;
  m_sortedBy = rhs.m_sortedBy;
  return *this;
}

#endif //DORQO
