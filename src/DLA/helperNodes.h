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
#include "tempVarNode.h"

class InputNode : public DLANode
{
  NodeType m_type;
 protected:
  DataTypeInfo m_dataTypeInfo;
#if TWOD
  const SizeList *m_msize;
  const SizeList *m_nsize;
#if DODM
  const SizeList *m_mlsize, *m_nlsize;
#endif
#else
  bool m_isScalar;
  SizesVec m_sizes;
#if DODM
  SizesVec m_lsizes;
#endif//DODM
#endif
  Name m_varName;

 public:
  InputNode();
#if TWOD
#if DOLLDLA
  Size m_rowStrideVal, m_colStrideVal;
  InputNode(NodeType type, Size m, Size n,
	    string name, 
	    Size rowStrideVal, Size colStrideVal,
	    string numRowsVar, string numColsVar,
	    string rowStrideVar, string colStrideVar,
	    Type dataType);

  explicit InputNode(string name, Size m, Size n,
	    Size rowStrideVal, Size colStrideVal,
	    Type dataType);

  string DataDeclaration();
  string RowStrideDefine();
  string ColStrideDefine();
  string NumRowsDefine();
  string NumColsDefine();

  bool IsContiguous();
#else
  InputNode(NodeType type, Size m, Size n, string name);
#endif
#if DODM
  InputNode(NodeType type, Size m, Size n, string name, DistType dist);
#endif
#elif DOTENSORS
  InputNode(NodeType type, const SizesVec &sizes, string name);
  InputNode(NodeType type, const SizesVec &sizes, const DistType &dist, string name);
  InputNode(NodeType type, const SizeList *sizes, string name);
#endif
  virtual NodeType GetType() const {return m_type;}
#if DOLLDLA
  virtual const Type GetDataType() const { return m_dataTypeInfo.m_type; }
#endif
  static Node* BlankInst() { return  new InputNode; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "inputNode";}
#if TWOD
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;
#if DODM
  virtual const SizeList* LocalM(ConnNum num) const;
  virtual const SizeList* LocalN(ConnNum num) const;
#endif
#elif DOTENSORS
  virtual const Dim NumDims(ConnNum num) const;
  virtual const SizeList* Len(ConnNum num, Dim dim) const;
  virtual const SizeList* LocalLen(ConnNum num, Dim dim) const;
#endif
  virtual Name GetName(ConnNum num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
};

class OutputNode : public DLANode
{
 public:
  virtual NodeType GetType() const {return "out";}
  static Node* BlankInst() { return  new OutputNode; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out) {}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "outputNode";}
#if TWOD
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;
#if DODM
  virtual const SizeList* LocalM(ConnNum num) const;
  virtual const SizeList* LocalN(ConnNum num) const;
#endif
#elif DOTENSORS
  virtual const Dim NumDims(ConnNum num) const;
  virtual const SizeList* Len(ConnNum num, Dim dim) const;
  virtual const SizeList* LocalLen(ConnNum num, Dim dim) const;
#endif
  virtual Name GetName(ConnNum num) const;
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}

#if DOTENSORS
  virtual bool CreatesNewVars() const {return false;}
#endif
};

#if DOLLDLA
//Constant value (e.g. used for Axpy)
class ConstVal : public DLANode
{
  Name m_varName;
  Coef m_val;
  const SizeList *m_sizes;
 public:
  ConstVal(string name, Coef val);
  virtual NodeType GetType() const {return "const val";}
  static Node* BlankInst() { return  new ConstVal("constVal", COEFZERO); }
  virtual const DataTypeInfo& DataType(ConnNum num) const {throw;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "const val";}
#if TWOD
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;
#if DODM
  virtual const SizeList* LocalM(ConnNum num) const {return ONES;}
  virtual const SizeList* LocalN(ConnNum num) const {return ONES;}
#endif
#else
blah
#endif
  virtual void ClearDataTypeCache();
 virtual void BuildDataTypeCache();
  virtual ~ConstVal();
 virtual Name GetName(ConnNum num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
};
#endif


#if TWOD
#if DOELEM
class MakeTrapNode : public DLAOp<1,1>
{
  Side m_side;
  Tri m_tri;
  int m_offset;
 public:
 MakeTrapNode(Side side, Tri tri, int offset) : m_side(side), m_tri(tri), m_offset(offset) {}
  static Node* BlankInst() { return  new MakeTrapNode(LEFT,LOWER,0); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const {return "MakeTrap";}
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "makeTrap";}
  virtual bool CanTrans() const;
  virtual bool CanTransposeInputs() const {return true;}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
};

class MoveMakeTrap : public SingleTrans
{
 public:
  bool CanApply(const Node *node) const;
  void Apply(Node *node) const;
};
#endif //DOELEM
#endif //TWOD

#if DOELEM||DOBLIS||DOTENSORS
class ScaleNode : public DLAOp<1,1>
{
 public:
  Coef m_val;
  ScaleNode(Layer layer, Coef val);
  static Node* BlankInst() { return  new ScaleNode(ABSLAYER, COEFZERO); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual NodeType GetType() const {return "Scale";}
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "scale";}
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void Prop();
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
};

class RemoveScaleByOne : public SingleTrans
{
 public:
  RemoveScaleByOne() {}
  virtual string GetType() const { return "RemoveScaleByOne"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return false;}
};
#endif //DOELEM||DOBLIS||DOTENSORS


#if TWOD
#if DOELEM||DOBLIS
class ScaleTrapNode : public DLAOp<1,1>
{
  Side m_side;
  Tri m_tri;
  Coef m_val;
 public:
 ScaleTrapNode(Layer layer, Side side, Tri tri, Coef val) : m_side(side), m_tri(tri), m_val(val) {SetLayer(layer);}
  static Node* BlankInst() { return  new ScaleTrapNode(ABSLAYER,LEFT,LOWER,COEFONE); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual NodeType GetType() const {return "ScaleTrap";}
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "scaleTrap";}
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void Prop();
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
};
#endif

#if DOELEM
//View a panel from two inputs
class ViewPan : public DLANode
{
 public:
  bool m_isVert;
  string m_name;
  SizeList *m_sizes;
#if DODM
  SizeList *m_lsizes;
#endif
 ViewPan(bool isVert, string name) 
   : m_isVert(isVert), m_name(name), m_sizes(NULL)
#if DODM
    , m_lsizes(NULL) 
#endif
    {}
  virtual ~ViewPan();
  virtual void Prop();
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual unsigned int NumOutputs() const {return 1;}
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
#if TWOD
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;
#if DODM
  virtual const SizeList* LocalM(ConnNum num) const;
  virtual const SizeList* LocalN(ConnNum num) const;
#endif
#else
sdlkfj
#endif
  virtual Name GetName(ConnNum num) const;
  virtual NodeType GetType() const {return "ViewPan";}
  static Node* BlankInst() { return new ViewPan(false,""); }
  virtual const DataTypeInfo& DataType(ConnNum num) const {return InputDataType(0);}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ViewPan";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
};


#if DOBLIS||DOELEM
//View the square around the diagonal on block
class ViewAroundDiag : public DLANode
{
 public:
  bool m_isVert;
  string m_name;
  Sizes *m_sizes0,*m_sizes1;
#if DODM 
  Sizes *m_lsizes0,*m_lsizes1;
#endif
 ViewAroundDiag(bool isVert, string name) 
   : m_isVert(isVert), m_name(name), 
    m_sizes0(NULL), m_sizes1(NULL)
#if DODM
    , m_lsizes0(NULL), m_lsizes1(NULL) 
#endif
    {}
  virtual ~ViewAroundDiag();
  virtual void Prop();
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual unsigned int NumOutputs() const {return 2;}
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
#if TWOD
  virtual const Sizes* GetM(ConnNum num) const;
  virtual const Sizes* GetN(ConnNum num) const;
#if DODM
  virtual const Sizes* LocalM(ConnNum num) const;
  virtual const Sizes* LocalN(ConnNum num) const;
#endif
#else
sdlkjf
#endif
  virtual Name GetName(ConnNum num) const;
  virtual NodeType GetType() const {return "ViewAround";}
  static Node* BlankInst() { return  new ViewAroundDiag(false,""); }
  virtual const DataTypeInfo& DataType(ConnNum num) const {return InputDataType(0);}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ViewAroundDiag";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
};
#endif //DOBLIS||DOELEM

#if DOBLIS||DOELEM
//Recombine the square around the diagonal and the other
// part of a block
class ViewAroundDiagCombine : public DLAOp<5,3>
{
 public:
  bool m_isVert;
 ViewAroundDiagCombine(bool isVert) : m_isVert(isVert) {}
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  static Node* BlankInst() { return  new ViewAroundDiagCombine(false); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void PrintCode(IndStream &out) {}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ViewAroundDiagCombine";}
  virtual NodeType GetType() const {return "ViewAroundCombine";}
  virtual void Prop();
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
};
#endif //DOBLIS||DOELEM
#endif


#if DOBLIS||DOELEM
class ViewTL : public DLANode
{
 public:
  ViewTL(Layer layer) {SetLayer(layer);}
#if TWOD
  virtual const Sizes* GetM(ConnNum num) const
  { return GetInputM(2); }
  virtual const Sizes* GetN(ConnNum num) const
  { return GetInputN(1); }
#if DODM
  virtual const Sizes* LocalM(ConnNum num) const
  { return InputLocalM(2); }
  virtual const Sizes* LocalN(ConnNum num) const
  { return InputLocalN(1); }
#endif
#else
bljsdf
#endif
  virtual Name GetName(ConnNum num) const;
  virtual void Prop();
  virtual unsigned int NumOutputs() const {return 1;}
  static Node* BlankInst() { return  new ViewTL(ABSLAYER); }
  virtual const DataTypeInfo& DataType(ConnNum num) const {return InputDataType(0);}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ViewTL";}
  virtual NodeType GetType() const {return "ViewTL";}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
};


class ViewTLCombine : public DLANode
{
 public:
  ViewTLCombine(Layer layer) {SetLayer(layer);}

  virtual const Sizes* GetM(ConnNum num) const
  { return GetInputM(1); }
  virtual const Sizes* GetN(ConnNum num) const
  { return GetInputN(1); }
#if DODM
  virtual const Sizes* LocalM(ConnNum num) const
  { return InputLocalM(1); }
  virtual const Sizes* LocalN(ConnNum num) const
  { return InputLocalN(1); }
#endif
  virtual Name GetName(ConnNum num) const {return GetInputName(1);}
  virtual void Prop();
  virtual unsigned int NumOutputs() const {return 1;}
  static Node* BlankInst() { return  new ViewTLCombine(ABSLAYER); }
  virtual const DataTypeInfo& DataType(ConnNum num) const {return InputDataType(1);}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void PrintCode(IndStream &out) {}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ViewTLCombine";}
  virtual NodeType GetType() const {return "ViewTLCombine";}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
};
#endif //DOBLIS||DOELEM
#endif
