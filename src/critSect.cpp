/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2013, The University of Texas and Bryan Marker

    DxTer is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DxTer is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.               

    You should have received a copy of the GNU General Public License
    along with DxTer.  If not, see <http://www.gnu.org/licenses/>.
*/

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
