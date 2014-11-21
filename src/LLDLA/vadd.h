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

class VAdd : public DLAOp<2, 1>
{
 public:
  VecType m_vecType;

  VAdd(VecType vecType, Layer layer);

  virtual void PrintCode(IndStream &out);
  virtual void Prop();
  virtual Phase MaxPhase() const;

  static Node* BlankInst();
  virtual Node* GetNewInst();

  static ClassType GetClass() { return "LLDLAVAdd"; }
  virtual ClassType GetNodeClass() const { return GetClass(); }

  virtual NodeType GetType() const;

 private:
  void PrintRowStride(IndStream &out);
  void PrintColStride(IndStream &out);
  void PrintGeneralStride(IndStream &out);
  void VectorOpInputDimensionCheck(ConnNum inputNum);
};

// BAMTODO: THIS IS PART OF A TEST
// TO SEE IF SPLIT / COMBINE ACTUALLY WORKS
/*
class SplitUnusualSizeVAdd : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  unsigned int m_multiple;

  SplitUnusualSizeVAdd(Layer fromLayer, Layer toLayer, unsigned int m_multiple);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const { return true; }
};
*/
class VAddLoopRef : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  VecType m_vtype;
  BSSize m_bs;

  VAddLoopRef(Layer fromLayer, Layer toLayer, VecType vtype, BSSize bs); 
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
