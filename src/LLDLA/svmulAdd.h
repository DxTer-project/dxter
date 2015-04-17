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

#include "LLDLA.h"

#if DOLLDLA

#include "DLAOp.h"

class SVMulAdd : public DLAOp<3, 1>
{
 private:
  VecType m_vecType;

 public:
  SVMulAdd(Layer layer, VecType vecType);

  virtual void PrintCode(IndStream &out) { throw; }
  virtual void Prop();
  virtual Phase MaxPhase() const;

  static Node* BlankInst();
  virtual Node* GetNewInst();
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);

  static ClassType GetClass() { return "LLDLASVMulAdd"; }
  virtual ClassType GetNodeClass() const { return GetClass(); }

  virtual NodeType GetType() const;
  virtual VecType GetVecType() const;

};

#endif // DOLLDLA
