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

#include "transform.h"
#include "DLAOp.h"

Loop* Chol1LowerAlg(Node *in, unsigned int num, bool dist);

Loop* Chol2LowerAlg(Node *in, unsigned int num, bool dist);
Loop* Chol2UpperAlg(Node *in, unsigned int num, bool dist);

Loop* Chol3LowerAlg(Node *in, unsigned int num, bool dist);
Loop* Chol3UpperAlg(Node *in, unsigned int num, bool dist);

class Chol : public DLAOp<1,1>
{
 public:
  Tri m_tri;
 Chol(Layer layer, Tri tri) : m_tri(tri) { SetLayer(layer); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void Prop();
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "chol";}
  virtual NodeType GetType() const {return "chol " + LayerNumToStr(m_layer);}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool ShouldCullDP() const;
  virtual Phase MaxPhase() const;
  static Node* BlankInst() { return new Chol(ABSLAYER, LOWER); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void SanityCheck();
  virtual void PrintCode(IndStream &out);
#if DOELEM
  virtual DistType GetDistType(unsigned int num) const;
#endif
};

class CholLoopExp : public SingleTrans
{
 public:
  int m_varNum;
  CholLoopExp(int varNum) : m_varNum(varNum) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

class DistCholToLocalChol : public SingleTrans
{
 public:
  virtual string GetType() const {return "Distributed Cholesky to Local Cholesky";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};
