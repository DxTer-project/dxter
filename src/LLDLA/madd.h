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

#include "DLAOp.h"
#include "LLDLA.h"
#include "loopSupport.h"

#if DOLLDLA

class MAdd : public DLAOp<2, 1>
{
 public:
  MAdd(Layer layer);
  virtual void PrintCode(IndStream &out);
  virtual void Prop();
  virtual Phase MaxPhase() const;

  static Node* BlankInst();
  virtual Node* GetNewInst() { return BlankInst(); }

  static ClassType GetClass() { return "LLDLAPrimMAdd"; }
  virtual ClassType GetNodeClass() const { return GetClass(); }

  virtual NodeType GetType() const;

 private:
  void PrintRowStride(IndStream &out);
  void PrintColStride(IndStream &out);
  void PrintGeneralStride(IndStream &out);

};

class MAddLoopRef : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  DimName m_dim;
  VecType m_vtype;
  BSSize m_bs;

  MAddLoopRef(Layer fromLayer, Layer toLayer, DimName dim, BSSize bs);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const { return true; }
};

class MAddToVAddLoopRef : SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  DimName m_dim;
  VecType m_vtype;
  BSSize m_bs;

  MAddToVAddLoopRef(Layer fromLayer, Layer toLayer, VecType vtype, BSSize bs);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const { return true; }
};

class MAddToRegArith : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  DimName m_dim;

  MAddToRegArith(Layer fromLayer, Layer toLayer);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const { return true; }
};

class MAddLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Size m_bs;

  MAddLowerLayer(Layer fromLayer, Layer toLayer, Size bs);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

#endif //DOLLDLA
