#include "critSect.h"
#include "loopSupport.h"

void CritSect::SanityCheck()
{
  if (!m_ownerPoss->m_pset->IsLoop())
    throw;
  PSet::SanityCheck();
}

void CritSect::PrintCurrPoss(IndStream &out, unsigned int &graphNum)
{
  if (!m_ownerPoss->m_pset->IsLoop())
    throw;
  Loop *loop = (Loop*)(m_ownerPoss->m_pset);
  Comm comm = loop->m_comm;
  if (comm == CORECOMM)
    throw;
  *out << "Critical section with communicator " << CommToStr(comm) << "; need correct output code\n";
  out.Indent();
  *out << "GetMutex(" << CommToStr(comm) << ");\n";

  PSet::PrintCurrPoss(out, graphNum);
  
  out.Indent();
  *out << "ReleaseMutex(" << CommToStr(comm) << ");\n";
}
