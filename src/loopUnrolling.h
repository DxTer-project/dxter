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

#include "layers.h"

#if DOLLDLA

#include "realLoop.h"
#include "shadowLoop.h"

class FullyUnrollLoop : public SingleTrans
{
 public:
  int m_numIters;
  FullyUnrollLoop(int numIters);

  virtual string GetType() const {return "fullyUnrollLoop " + std::to_string((long long int) m_numIters);}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class ViewMultipleIters : public DLANode
{
 public:
  PartDir m_partDir;
  BSSize m_bs;
  int m_numIters;
  Sizes *m_sizes;
  DataTypeInfo m_info;
 ViewMultipleIters(PartDir partDir, BSSize bs, int numIters)
   : m_partDir(partDir), m_bs(bs), m_numIters(numIters),
    m_sizes(NULL)
    {}
  virtual ~ViewMultipleIters();
  virtual void Prop();
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual unsigned int NumOutputs() const {return m_numIters+1;}
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
  virtual const Sizes* GetM(ConnNum num) const;
  virtual const Sizes* GetN(ConnNum num) const;
  virtual Name GetName(ConnNum num) const;
  virtual void AddVariables(VarSet &set) const;
  virtual NodeType GetType() const {return "ViewMultipleIters";}
  static Node* BlankInst() { return  new ViewMultipleIters(LASTPARTDIR, (BSSize)BADBSSIZE, -1); }
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  bool KeepsInputVarLive(Node *input, ConnNum numIn, ConnNum &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ViewMultipleIters";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
};

class CombineMultipleIters : public DLANode
{
 public:
  PartDir m_partDir;
  BSSize m_bs;
  int m_numIters;
 CombineMultipleIters(PartDir partDir, BSSize bs, int numIters)
   : m_partDir(partDir), m_bs(bs), m_numIters(numIters)
    {}
  virtual void Prop();
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual unsigned int NumOutputs() const {return 1;}
  virtual const Sizes* GetM(ConnNum num) const;
  virtual const Sizes* GetN(ConnNum num) const;
  virtual Name GetName(ConnNum num) const;
  virtual NodeType GetType() const {return "CombineMultipleIters";}
  static Node* BlankInst() { return  new CombineMultipleIters(LASTPARTDIR, (BSSize)BADBSSIZE, -1); }
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  bool KeepsInputVarLive(Node *input, ConnNum numIn, ConnNum &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "CombineMultipleIters";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
};


#endif //DOLLDLA
