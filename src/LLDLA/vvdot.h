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

#include "DLAOp.h"
#include "loopSupport.h"
#include "LLDLA.h"

#if DOLLDLA

class VVDot : public DLAOp<3, 1>
{
 public:
  Type m_type;
  
  VVDot(Layer layer, Type type);

  virtual void PrintCode(IndStream &out);
  virtual void Prop();
  virtual Phase MaxPhase() const;

  static Node* BlankInst();
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);

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
 VVDotLoopRef(Layer fromLayer, Layer toLayer, BSSize bs) 
   : m_fromLayer(fromLayer), m_toLayer(toLayer), m_bs(bs) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const { return true; }
};

class VVDotLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Size m_bs;
 VVDotLowerLayer(Layer fromLayer, Layer toLayer, Size bs)
   :m_fromLayer(fromLayer), m_toLayer(toLayer), m_bs(bs) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

#endif // DOLLDLA
