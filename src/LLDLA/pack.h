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
#include "LLDLA.h"

#ifndef PACK_H_
#define PACK_H_

#if DOLLDLA

class Pack : public DLAOp<2, 1> {
 public:
  explicit Pack(Layer layer);

  virtual NodeType GetType() const { return "Pack"; }
  virtual ClassType GetNodeClass() const { return GetClass(); }
  static ClassType GetClass() { return "Pack"; }

  virtual const DataTypeInfo& DataType(ConnNum num) const;
  virtual bool Overwrites(const Node* input, ConnNum num) const;

  virtual bool IsReadOnly() const { return false; }
  virtual bool IsDataDependencyOfInput() const { return true; }

  /*  virtual Sizes* GetInputM(ConnNum num) const;
      virtual Sizes* GetInputN(ConnNum num) const;*/

  virtual int PackM();
  virtual int PackN();

  virtual Dir PackDir() = 0;

  virtual void Prop();
  virtual void PrintCode(IndStream& out);
};

#endif // DOLLDLA

#endif // PACK_H_
