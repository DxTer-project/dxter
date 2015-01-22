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

#include "DLAOp.h"
#include "loopSupport.h"
#include "LLDLA.h"

#if DOLLDLA

class VVDot : public DLAOp<3, 1>
{
 public:
  VVDot(Layer layer);

  virtual void PrintCode(IndStream &out);
  virtual void Prop();
  virtual Phase MaxPhase() const;

  static Node* BlankInst();
  virtual Node* GetNewInst() { return BlankInst(); }

  static ClassType GetClass() { return "LLDLASVVDot"; }
  virtual ClassType GetNodeClass() const { return GetClass(); }

  virtual NodeType GetType() const;

 private:
  void PrintRowStride(IndStream &out);
  void PrintColStride(IndStream &out);
  void PrintGeneralStride(IndStream &out);
  void VectorOpInputDimensionCheck(ConnNum inputNum);

};

class VVDotLoopRef : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  BSSize m_bs;
  VVDotLoopRef(Layer fromLayer, Layer toLayer, BSSize bs);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const { return true; }
};

class VVDotToRegArith : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  VVDotToRegArith(Layer fromLayer, Layer toLayer);
  virtual string GetType() const;
  virtual bool CanApply(const Node* node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const { return true; }
};

class VVDotLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Size m_bs;
  VVDotLowerLayer(Layer fromLayer, Layer toLayer, Size bs);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

#endif // DOLLDLA
