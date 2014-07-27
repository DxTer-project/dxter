/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2014, The University of Texas and Bryan Marker

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

#include "splitSingleIter.h"
#include "combineSingleIter.h"
#include "elemRedist.h"
#include <cmath>
#include "LLDLA.h"

SplitBase::SplitBase() : LoopTunnel(LASTTUNNEL), 
#if TWOD
m_dir(LASTPARTDIR)
#else
m_partDim(99)
#endif
{
}

#if TWOD
SplitBase::SplitBase(PartDir dir, TunType type, bool isControl) 
  : LoopTunnel(type), m_dir(dir), m_isControlTun(isControl)
{
}
#else
SplitBase::SplitBase(unsigned int partDim, TunType type, bool isControl) 
  : LoopTunnel(type), m_partDim(partDim), m_isControlTun(isControl)
{
}
#endif



void SplitBase::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  LoopTunnel::Duplicate(orig, shallow, possMerging);
  const SplitBase *split = (SplitBase*)orig;
#if TWOD
  m_dir = split->m_dir;
#else
  m_partDim = split->m_partDim;
#endif
  m_isControlTun = split->m_isControlTun;
}


void SplitBase::FlattenCore(ofstream &out) const
{
  LoopTunnel::FlattenCore(out);
#if TWOD
  WRITE(m_dir);
#else
  //no m_lsizes - in loop tunnel
  WRITE(m_partDim);
#endif
  WRITE(m_isControlTun);

}


void SplitBase::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  LoopTunnel::UnflattenCore(in,info);
#if TWOD
  READ(m_dir);
#else
  //no m_lsizes - in loop tunnel
  READ(m_partDim);
#endif
  READ(m_isControlTun);

}
