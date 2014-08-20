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

#if DOLLDLA

#include "DLAOp.h"
#include "LLDLA.h"
#include "loopSupport.h"
#include "regArith.h"
#include "regLoadStore.h"

class VMMul : public DLAOp<3, 1>
{
 public:
  Type m_type;
  int m_regWidth;
  VMMul(Layer layer, Type type);

  virtual void PrintCode(IndStream &out);
  virtual void Prop();
  virtual Phase MaxPhase() const;

  static Node* BlankInst();
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);

  static ClassType GetClass() { return "LLDLAPrimVMMul"; }
  virtual ClassType GetNodeClass() const { return GetClass(); }

  virtual NodeType GetType() const;

 private:
  void PrintRowStride(IndStream &out);
  void PrintColStride(IndStream &out);
  void PrintGeneralStride(IndStream &out);

};

class VMMulLoopRef : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  DimName m_dim;
  BSSize m_bs;
  Type m_type;
  int m_regWidth;

  VMMulLoopRef(Layer fromLayer, Layer toLayer, DimName dim, BSSize bs, Type type);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const { return true; }

 private:
  void ApplyDimN(Node* node) const;
  void ApplyDimK(Node* node) const;
};

class VMMulToRegArith : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Type m_type;
  int m_regWidth;
  VMMulToRegArith(Layer fromLayer, Layer toLayer, Type type);
  virtual string GetType() const;
  virtual bool CanApply(const Node* node) const;
  virtual void Apply(Node* node) const;
  virtual bool IsRef() const { return true; }
};

class VMMulLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Size m_bs;
  Type m_type;
  int m_regWidth;
  VMMulLowerLayer(Layer fromLayer, Layer toLayer, Size bs, Type type);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const { return true; }
};

#endif // DOLLDLA
