/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2015, The University of Texas and Bryan Marker

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



#include "combineSingleIter.h"
#include "splitSingleIter.h"
#include "elemRedist.h"
#include <cmath>

void CombineBase::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  LoopTunnel::Duplicate(orig, shallow, possMerging);
  const CombineBase *com = (CombineBase*)orig;
#if TWOD
  m_dir = com->m_dir;
#else
  m_partDim = com->m_partDim;
#endif
}
void CombineBase::FlattenCore(ofstream &out) const
{
  LoopTunnel::FlattenCore(out);
#if TWOD
  WRITE(m_dir);
#else
  WRITE(m_partDim);
#endif
}

void CombineBase::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  LoopTunnel::UnflattenCore(in,info);
#if TWOD
  READ(m_dir);
#else
  READ(m_partDim);
#endif
}
