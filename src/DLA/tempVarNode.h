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

#include "DLANode.h"
#include "DLAOp.h"
#include "LLDLA.h"

class TempVarNode : public DLANode
{
 public:
  DataTypeInfo m_info;
  string m_name;
#if DOELEM
  Sizes *m_mlsize, *m_nlsize;
#elif DOTENSORS
  SizesArray m_lsizes;
  SizesArray m_sumLens;
  EntryList m_sumDims;
  Sizes m_ones;
#endif

  TempVarNode();
  TempVarNode(string name);
#if DOTENSORS
  TempVarNode(DistType dist, EntryList sumDims);
  TempVarNode(DistType dist, EntryList sumDims, string name);
#endif
#if DODM
  TempVarNode(DistType dist);
  TempVarNode(DistType dist, string name); 
#endif

 ~TempVarNode();
  virtual NodeType GetType() const;
  static Node* BlankInst() { return  new TempVarNode; }
  bool KeepsInputVarLive(Node *input, ConnNum numIn, ConnNum &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);

  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "tempVar";}
  virtual const DataTypeInfo& DataType(ConnNum num) const;
#if TWOD
  virtual const Sizes* GetM(ConnNum num) const;
  virtual const Sizes* GetN(ConnNum num) const;
#if DODM
  virtual const Sizes* LocalM(ConnNum num) const;
  virtual const Sizes* LocalN(ConnNum num) const;
#endif
#else
  virtual const Dim NumDims(ConnNum num) const;
  virtual const Sizes* Len(ConnNum num, Dim dim) const;
  virtual const Sizes* LocalLen(ConnNum num, Dim dim) const;
#endif
  virtual Name GetName(ConnNum num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
  virtual bool IsReadOnly() const {return true;}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual bool IsDataDependencyOfInput() const {return false;}
};

class MoveTempVarNodeIntoLoop : public SingleTrans
{
 public:
  virtual string GetType() const { return "MoveTempVarNodeIntoLoop"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};
