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

#include "transform.h"
#include "DLAOp.h"
#include "elemRedist.h"
#include "gemm.h"

#if DOELEM

RealLoop* HetrmmAlgVar1Lower(Node *in, ConnNum num);
RealLoop* HetrmmAlgVar1Upper(Node *in, ConnNum num);

RealLoop* HetrmmAlgVar2Lower(Node *in, ConnNum num);
RealLoop* HetrmmAlgVar2Upper(Node *in, ConnNum num);

RealLoop* HetrmmAlgVar3Lower(Node *in, ConnNum num);
RealLoop* HetrmmAlgVar3Upper(Node *in, ConnNum num);


class Hetrmm : public DLAOp<1,1>
{
 public:
  Tri m_tri;
 Hetrmm(Layer layer, Tri tri) : m_tri(tri) {SetLayer(layer);}
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Hetrmm";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual NodeType GetType() const;
  static Node* BlankInst() { return new Hetrmm(ABSLAYER, LOWER); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual Phase MaxPhase() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual bool ShouldCullDP() const;
};

class HetrmmLoopExp : public SingleTrans
{
 public:
  int m_var;
 HetrmmLoopExp(int var) : m_var(var) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

#if DODPPHASE
class DistHetrmmToLocalHetrmm : public SingleTrans
{
 public:
  virtual string GetType() const {return "Distributed Hetrmm to Local Hetrmm";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};
#endif

#endif //DOBLIS||DOELEM
