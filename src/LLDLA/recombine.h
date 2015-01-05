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
#include "node.h"
#include "layers.h"

#if DOLLDLA

class Recombine : public DLAOp<3, 1>
{

 private:
  Dir m_partType;

 public:
  Layer m_layer;
  
  Recombine(Layer layer, Dir partType);

  virtual void PrintCode(IndStream &out);

  virtual ~Recombine() {}
  inline void SetLayer(Layer layer) { m_layer = layer; }
  inline Layer GetLayer() const { return m_layer; }
  virtual bool IsReadOnly() const { return false; }
  virtual bool CanTrans() const { return false; }

  virtual NodeType GetType() const { return "recombineNode"; }
  static ClassType GetClass() { return "recombineNode"; }
  virtual ClassType GetNodeClass() const { return GetClass(); }

  virtual Node* GetNewInst() { return BlankInst(); }
  static Node* BlankInst() { return new Recombine(ABSLAYER, VERTICAL); }
  virtual void Duplicate(const Node* orig, bool shallow, bool possMerging);

  virtual Name GetName(ConnNum num) const;

  virtual void Prop();

  virtual const Sizes* GetM(ConnNum num) const;
  virtual const Sizes* GetN(ConnNum num) const;

  virtual const DataTypeInfo& DataType(ConnNum num) const;
};

class RecombineLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Size m_bs;

  RecombineLowerLayer(Layer fromLayer, Layer toLayer);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

#endif //DOLLDLA
