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

#include "layers.h"

#if DOLLDLA

#include "gemmTransformations.h"
#include "realLoop.h"

class MMulLoopExp : public GemmLoopExp
{
 protected:
  bool CheckGemmLoop(const Node* node) const;

 public:
  Type m_type;
  MMulLoopExp(Layer fromLayer, Layer toLayer, DimName dim, BSSize bsSize, Type type);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class GemmTransToNotTrans : public SingleTrans
{
 public:
  Layer m_layer;
  Type m_type;
 GemmTransToNotTrans(Layer layer, Type type) :
  m_layer(layer), m_type(type) {}
  virtual string GetType() const {return "GemmTransToNotTrans";}
  virtual bool IsRef() const {return true;}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class MMulToMVMul : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Size m_bs;
  Type m_type;
 MMulToMVMul(Layer fromLayer, Layer toLayer)
   : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}  
};

class MMulToVMMul : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Size m_bs;
  Type m_type;
 MMulToVMMul(Layer fromLayer, Layer toLayer)
   : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}  
};

class LLDAGemmLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Size m_bs;
  Type m_type;
 LLDAGemmLowerLayer(Layer fromLayer, Layer toLayer, Size bs, Type type)
   : m_fromLayer(fromLayer), m_toLayer(toLayer), m_bs(bs), m_type(type) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

#endif
