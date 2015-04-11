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

#include "DLANode.h"
#include "DLAOp.h"
#include "LLDLA.h"
#include "sizes.h"

class TempVarNode : public DLANode
{
 public:
  DataTypeInfo m_info;
  string m_name;
#if DOELEM
  Sizes *m_mlsize, *m_nlsize;
#elif DOTENSORS
  SizesVec m_lsizes;
  SizesVec m_sumLens;
  EntryList m_sumDims;
  const SizeList *m_ones;
  string m_align;
  DimVec m_alignModes;
  DimVec m_alignModesSrc;
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

  virtual NodeType GetType() const;
  static Node* BlankInst() { return  new TempVarNode; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);

  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "tempVar";}
  virtual const DataTypeInfo& DataType(ConnNum num) const;
#if TWOD
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;
#if DODM
  virtual const SizeList* LocalM(ConnNum num) const;
  virtual const SizeList* LocalN(ConnNum num) const;
#endif
#else
  virtual const Dim NumDims(ConnNum num) const;
  virtual const SizeList* Len(ConnNum num, Dim dim) const;
  virtual const SizeList* LocalLen(ConnNum num, Dim dim) const;
#endif
  virtual Name GetName(ConnNum num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
  virtual bool IsReadOnly() const {return true;}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual bool IsDataDependencyOfInput() const {return false;}
#if DOTENSORS
  virtual void AddVariables(VarSet &set) const;
#endif
};

class MoveTempVarNodeIntoLoop : public SingleTrans
{
 public:
  virtual string GetType() const { return "MoveTempVarNodeIntoLoop"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};


class MoveTempVarNodeIntoSet : public SingleTrans
{
 public:
  virtual string GetType() const { return "MoveTempVarNodeIntoSet"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class TempVarFromTempVar : public SingleTrans
{
 public:
  virtual string GetType() const { return "TempVarFromTempVar"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};
