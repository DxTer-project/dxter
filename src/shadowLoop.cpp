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



#include "loopSupport.h"
#include "elemRedist.h"
#include <cmath>
#include <climits>
#include "pack.h"
#include "critSect.h"
#include "blis.h"

int ShadowLoop::GetBS() const
{
  return GetBSSize().GetSize();
}

BSSize ShadowLoop::GetBSSize() const
{
  return ((RealLoop*)m_realPSet)->GetBSSize();
}

#if TWOD
DimName ShadowLoop::GetDimName() const
{
  return ((RealLoop*)m_realPSet)->GetDimName();
}
#endif

void ShadowLoop::Prop()
{
  if (!BasePSet::m_hasProped) {
    if (!m_realPSet || !m_realPSet->IsLoop())
      throw;
    IntLoop<ShadowPSet>::Prop();
  }
}

LoopType ShadowLoop::GetType() const
{
  return ((RealLoop*)m_realPSet)->GetType();
}

ShadowPSet* ShadowLoop::GetNewShadow()
{
  return ((RealLoop*)m_realPSet)->GetNewShadow();
}


#if DOBLIS
bool ShadowLoop::HasIndepIters() const
{
  return m_realPSet->HasIndepIters();
}

bool IsParallel() const 
{
  return m_realPSet->IsParallel();
}
#endif



const IntSet& ShadowLoop::GetLabel() const
{
  return ((RealLoop*)m_realPSet)->GetLabel();
}
