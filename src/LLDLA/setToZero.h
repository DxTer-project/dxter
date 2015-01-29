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

#include "DLAOp.h"

#if DOLLDLA

class SetToZero : public DLAOp<1, 1>
{
 private:
  bool InputIsGenStride();

  void NaiveZeroPrintout(IndStream& out);
  void MemsetPrintout(IndStream& out);

 public:
  SetToZero(Layer layer);
  virtual NodeType GetType() const;
  static Node* BlankInst();
  bool KeepsInputVarLive(Node* input, ConnNum numIn, ConnNum& numOut) const;
  virtual Node* GetNewInst();

  virtual Phase MaxPhase() const;

  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const;
  static ClassType GetClass();

  virtual Name GetName(ConnNum num) const;

  virtual bool IsReadOnly() const;
  virtual bool Overwrites(const Node* input, ConnNum num) const;
  virtual bool IsDataDependencyOfInput() const;

};

class SetToZeroLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;

  SetToZeroLowerLayer(Layer fromLayer, Layer toLayer);

  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const;
};

#endif // DOLLDLA
