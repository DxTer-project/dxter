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

#ifndef _VADD_H_
#define _VADD_H_

#include "DLAOp.h"
#include "LLDLA.h"
#include "loopSupport.h"

#if DOLLDLA

string VecTypeToString(VecType vType);

class VAdd : public DLAOp<2, 1>
{
 private:
  VecType m_vecType;

 public:
  VAdd(Layer layer, VecType vecType);

  virtual void PrintCode(IndStream &out);
  virtual void Prop();
  virtual Phase MaxPhase() const;

  static Node* BlankInst();
  virtual Node* GetNewInst();
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);

  static ClassType GetClass() { return "LLDLAVAdd"; }
  virtual ClassType GetNodeClass() const { return GetClass(); }

  virtual NodeType GetType() const;
  virtual VecType GetVecType() const;

 private:
  void PrintRowStride(IndStream &out);
  void PrintColStride(IndStream &out);
  void PrintGeneralStride(IndStream &out);
  void VectorOpInputDimensionCheck(ConnNum inputNum);
};

class VAddLoopRef : public SingleTrans
{
 private:
  bool CheckRowVectorDimension(const VAdd* vadd) const;
  bool CheckColVectorDimension(const VAdd* vadd) const;

 public:
  Layer m_fromLayer, m_toLayer;
  VecType m_vtype;

  VAddLoopRef(Layer fromLayer, Layer toLayer, VecType vtype);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const { return true; }
};


class VAddLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Size m_bs;

  VAddLowerLayer(Layer fromLayer, Layer toLayer, Size bs);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

class VAddToRegArith : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  DimName m_dim;

  VAddToRegArith(Layer fromLayer, Layer toLayer);
  virtual string GetType() const;
  virtual bool CanApply(const Node* node) const;
  virtual void Apply(Node* node) const;
  virtual bool IsRef() const { return true; }
};

#endif // DOLLDLA

#endif
