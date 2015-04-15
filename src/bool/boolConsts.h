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
#if DOBOOL

#include "DLAOp.h"

class True : public Node
{
  static DataTypeInfo m_info;
  static Name m_name;
 public:
  static Node* BlankInst() { return new True; }
  virtual Node* GetNewInst() { return BlankInst(); }
  //  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "True";}
  virtual Cost GetCost() {return 0;}
  virtual Name GetName(ConnNum num) const {return m_name;}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual const DataTypeInfo& DataType(ConnNum num) const { return m_info;}
  virtual NodeType GetType() const {return "True";}
  virtual Phase MaxPhase() const {return NUMPHASES;}
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  //  virtual void AddVariables(VarSet &set) const;
};


class False : public Node
{
  static DataTypeInfo m_info;
  static Name m_name;
 public:
  static Node* BlankInst() { return new False; }
  virtual Node* GetNewInst() { return BlankInst(); }
  //  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "False";}
  virtual Cost GetCost() {return 0;}
  virtual Name GetName(ConnNum num) const {return m_name;}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual const DataTypeInfo& DataType(ConnNum num) const { return m_info;}
  virtual NodeType GetType() const {return "False";}
  virtual Phase MaxPhase() const {return NUMPHASES;}
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  //  virtual void AddVariables(VarSet &set) const;
};


#endif //DOBOOL
