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



#pragma once

#include "loop.h"
#include "loopTunnel.h"

class Combine : public LoopTunnel
{
 public:
#if TWOD
  PartDir m_dir;
#else
  Dim m_partDim;
#endif
#if TWOD
 Combine() :LoopTunnel(LASTTUNNEL),m_dir(LASTPARTDIR) {}
 Combine(PartDir dir, PossTunType type) : LoopTunnel(type), m_dir(dir) {}
#else
 Combine() :LoopTunnel(LASTTUNNEL),m_partDim(99) {}
   Combine(Dim partDim, PossTunType type) : LoopTunnel(type), m_partDim(partDim) {}
#endif
  virtual void PrintCode(IndStream &out);
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  static Node* BlankInst() { return new Combine(LASTPARTDIR,LASTTUNNEL);}
  virtual Node* GetNewInst() {return BlankInst(); }
  virtual NodeType GetType() const;
  virtual void SanityCheck();
  virtual void Prop();
  virtual PossTunnel* GetSetTunnel();
  virtual unsigned int NumOutputs() const {return 1;}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "combine";}
#if TWOD
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
#if DOLLDLA
  virtual Stride RowStride(unsigned int num) const;
  virtual Stride ColStride(unsigned int num) const;
#endif //DOLLDLA
#else
  virtual const Dim NumDims(unsigned int num) const;
  virtual const Sizes* Len(unsigned int num, Dim dim) const;
  virtual const Sizes* LocalLen(unsigned int num, Dim dim) const;
#endif
  virtual Name GetName(unsigned int num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
};
