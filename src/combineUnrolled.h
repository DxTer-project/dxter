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
#include "combineBase.h"

class CombineUnrolled : public CombineBase
{
 public:
  unsigned int m_unrollFactor;
#if TWOD
  CombineUnrolled();
  CombineUnrolled(PartDir dir, unsigned int unrollFactor, TunType type);
#else
  CombineUnrolled();
  CombineUnrolled(Dim partDim, unsigned int unrollFactor, TunType type);
#endif
  virtual void PrintCode(IndStream &out);
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  static Node* BlankInst() { return new CombineUnrolled(LASTPARTDIR,0,LASTTUNNEL);}
  virtual Node* GetNewInst() {return BlankInst(); }
  virtual NodeType GetType() const;
  virtual void Prop();
  virtual Tunnel* GetSetTunnel();
  virtual unsigned int NumOutputs() const {return 1;}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "combineUnrolled";}
  virtual const DataTypeInfo& DataType(ConnNum num) const;
#if TWOD
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;
#if DODM
  virtual const SizeList* LocalM(ConnNum num) const;
  virtual const SizeList* LocalN(ConnNum num) const;
#endif
#else
  virtual const Dim NumDims(ConnNum num) const;
  virtual const SizeList* Len(ConnNum num, Dim dim) const;
  virtual const SizeList* LocalLen(ConnNum num, Dim dim) const;
#endif
  virtual Name GetName(ConnNum num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool IsCombine() const {LOG_FAIL("replacement for throw call");return true;}
};
