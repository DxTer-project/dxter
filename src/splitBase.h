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

//LoopTunnel for spliting/indxing into a matrix
class SplitBase : public LoopTunnel
{
 public:
#if TWOD
  PartDir m_dir;
#else
  //no m_lsizes - in loop tunnel
  Dim m_partDim;
#endif
  bool m_isControlTun;
  SplitBase();
#if TWOD
  SplitBase(PartDir dir, TunType type, bool isControl);
#else
  SplitBase(unsigned int partDim, TunType type, bool isControl);
#endif
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool IsSplit() const {return true;}
  virtual unsigned int NumberOfLoopExecs() const = 0;
  virtual void PrintVarDeclarations(BSSize bs, IndStream &out) const = 0;
  virtual void PrintIncrementAtEndOfLoop(BSSize bs, IndStream &out) const = 0;
#if TWOD
  virtual unsigned int NumIters(Size bs, Size m, Size n) const = 0;
#else
  virtual unsigned int NumIters(Size bs, Size size) const = 0;
#endif
  virtual unsigned int NumIters(unsigned int iterNum) const = 0;
#if DOLLDLA
  virtual string LoopBound() = 0;
#endif
};
