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



#pragma once

#include "realLoop.h"
#include "shadowLoop.h"
#include "loopTunnel.h"

#if DOLOOPS
class CombineBase : public LoopTunnel
{
 public:
#if TWOD
  PartDir m_dir;
#elif DOTENSORS
  Dim m_partDim;
#endif
#if TWOD
 CombineBase() :LoopTunnel(LASTTUNNEL),m_dir(LASTPARTDIR) {}
 CombineBase(PartDir dir, TunType type) : LoopTunnel(type), m_dir(dir) {}
#elif DOTENSORS
 CombineBase() :LoopTunnel(LASTTUNNEL),m_partDim(99) {}
   CombineBase(Dim partDim, TunType type) : LoopTunnel(type), m_partDim(partDim) {}
#endif
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool IsCombine() const {LOG_FAIL("replacement for throw call"); return true;}
};
#endif //DOLOOPS
