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


#pragma once

#include "layers.h"

#if DOTENSORS

#include "DLANode.h"
#include "transform.h"
#include "DLAOp.h"

class SumScatterUpdateNode : public DLAOp<2,1>
{
 public:
  Coef m_coef;
  EntrySet m_sumDims;
 SumScatterUpdateNode() : DLAOp<2,1>(), m_coef(COEFVALZERO) {}
  SumScatterUpdateNode(Coef coeff, const EntrySet &sumDims);
  static Node* BlankInst() { return  new SumScatterUpdateNode; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual bool IsRedistNode() const {return true;}
  virtual NodeType GetType() const;
  virtual Phase MaxPhase() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  Dim MySumDim() const;
  //  RedistNode* CreateTrans(Trans trans);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "sumScatter";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual const DistType& GetDistType(unsigned int num) const { return InputDistType(1); }
  void CheckSumDimsInOutput() const;
};

class SeparateRedistFromSumScatter : public SingleTrans
{
 public:
  virtual string GetType() const { return "SeparateRedistFromSumScatter"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

class SplitSumScatter : public SingleTrans
{
 public:
  Dim m_dim;
 SplitSumScatter(Dim dim) : m_dim(dim) {}
  virtual string GetType() const { return (string)"SplitSumScatter" + (char)(m_dim+48); }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};


class MoveSumScatterRedistAfter : public SingleTrans
{
 public:
  virtual string GetType() const { return "MoveSumScatterRedistAfter"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};


#endif
