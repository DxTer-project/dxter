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

/*
  This file defines the LLDLAGemm class. This class is used in
  LLDLA operations to represent multiplication a gemm  whose input's
  row and column dimensions are both equal to mu.

  The only significant differences from Gemm are that this class appears
  only in the primitive layer of LLDLA operations and that it has 2
  fields (m_rowStride, m_colStride) which denotes the row and column
  stride of its arguments. Based on these two field's
  values LLDLAGemm's PrintCode method generates the appropriate API
  call.
*/

#pragma once

#include "gemm.h"

#if DOLLDLA

class LLDLAGemm : public Gemm
{
 public:
  int m_regWidth;
  LLDLAGemm(Coef alpha, Coef beta, Type type, Layer layer);
  virtual void PrintCode(IndStream &out);
  virtual void Prop();
  virtual Phase MaxPhase() const { return NUMPHASES; }

  static Node* BlankInst();
  virtual Node* GetNewInst() { return BlankInst(); }

  static ClassType GetClass() { return "LLDLAPrimGemm"; }
  virtual ClassType GetNodeClass() const { return GetClass(); }

  virtual NodeType GetType() const;

  void PrintRowStride(IndStream &out);
  void PrintColStride(IndStream &out);
  void PrintGeneralStride(IndStream &out);
};

class LLDLAGemmToPrim : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Type m_type;
  int m_regWidth;
  LLDLAGemmToPrim(Layer fromLayer, Layer toLayer, Type type);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};



#endif //DOLLDLA
