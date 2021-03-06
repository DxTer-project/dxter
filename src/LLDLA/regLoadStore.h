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

#if DOLLDLA

#include "loopTunnel.h"

class LoadToRegs : public DLANode
{
 protected:
  void SetCost();
  void SanityCheckInputDims();

 public:

  virtual NodeType GetType() const { return "LoadToRegs"; }
  static Node* BlankInst() { return  new LoadToRegs(); }
  virtual Node* GetNewInst() { return BlankInst(); }

  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "LoadtoRegs";}
  virtual const DataTypeInfo& DataType(ConnNum num) const {return InputDataType(0);}
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;

  virtual Name GetName(ConnNum num) const;

  virtual bool IsReadOnly() const {return true;}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual bool IsDataDependencyOfInput() const {return true;}

  virtual void AddVariables(VarSet &set) const;

  //  virtual Phase MaxPhase() const;
};

class PackedLoadToRegs : public DLANode {
 private:
  int ComputeResidual();
 public:

  virtual NodeType GetType() const { return "PackedLoadToRegs"; }
  static Node* BlankInst() { return  new PackedLoadToRegs(); }
  virtual Node* GetNewInst() { return BlankInst(); }

  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "PackedLoadtoRegs";}
  virtual const DataTypeInfo& DataType(ConnNum num) const {return InputDataType(0);}

  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual bool IsDataDependencyOfInput() const {return true;}
  virtual bool IsReadOnly() const {return true;}

  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;

  virtual Name GetName(ConnNum num) const;

  virtual void AddVariables(VarSet &set) const;
};

class DuplicateRegLoad : public DLANode
{
 public:
  const SizeList *m_mSizes;
  const SizeList *m_nSizes;
  DataTypeInfo m_info;

  virtual NodeType GetType() const { return "DuplicateRegLoad"; }
  static Node* BlankInst() { return  new DuplicateRegLoad(); }
  virtual Node* GetNewInst() { return BlankInst(); }

  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "DuplicateVecReg";}
  virtual const DataTypeInfo& DataType(ConnNum num) const {return m_info;}
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;

  virtual Name GetName(ConnNum num) const;

  virtual bool IsReadOnly() const {return true;}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual bool IsDataDependencyOfInput() const {return true;}

  virtual void AddVariables(VarSet &set) const;
  virtual void BuildDataTypeCache();
  virtual void ClearDataTypeCache();
};

class TempVecReg : public DLANode
{
 public:
  const SizeList *m_mSizes;
  const SizeList *m_nSizes;

  DataTypeInfo m_info;
  virtual NodeType GetType() const { return "TempVecReg"; }
  static Node* BlankInst() { return  new TempVecReg(); }
  virtual Node* GetNewInst() { return BlankInst(); }

  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "TempVecReg";}
  virtual const DataTypeInfo& DataType(ConnNum num) const {return m_info;}
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;

  virtual Name GetName(ConnNum num) const;

  virtual bool IsReadOnly() const {return true;}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual bool IsDataDependencyOfInput() const {return false;}

  virtual void AddVariables(VarSet &set) const;
  virtual void BuildDataTypeCache();
  virtual void ClearDataTypeCache();
};

class StoreFromRegs : public DLAOp<2,1>
{
 public:

  virtual NodeType GetType() const {return "StoreFromRegs";}
  static Node* BlankInst() { return  new StoreFromRegs(); }
  virtual Node* GetNewInst() { return BlankInst(); }

  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "StoreFromRegs";}

  virtual bool IsReadOnly() const {return false;}
  virtual bool IsDataDependencyOfInput() const {return true;}

  void StoreNonContigLocations(IndStream &out, string regVarName, string storePtr, string strideVar);
};

class UnpackStoreFromRegs : public DLAOp<2,1>
{
 private:
  int ComputeResidual();

 public:
  virtual NodeType GetType() const {return "UnpackStoreFromRegs";}
  static Node* BlankInst() { return  new UnpackStoreFromRegs(); }
  virtual Node* GetNewInst() { return BlankInst(); }

  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "UnpackStoreFromRegs";}

  virtual bool IsReadOnly() const {return false;}
  virtual bool IsDataDependencyOfInput() const {return true;}
};

class HoistLoad : public SingleTrans
{
 protected:
  bool OnlyFoundLoad(const LoopTunnel* setTunIn) const;
  bool PossChildrenAreOnlyLoads(const LoopTunnel* possTunIn) const;
  void RemoveLoadsFromPosses(LoopTunnel* setTunIn) const;
  void AddOneOuterLoad(LoopTunnel* setTunIn) const;

 public:
  virtual string GetType() const {return "HoistLoad";}
  virtual bool IsRef() const {return false;}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};


#endif //DOLLDLA
