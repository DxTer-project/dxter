/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2013, The University of Texas and Bryan Marker

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

class InputNode : public DLANode
{
  NodeType m_type;
#if TWOD
  Sizes m_msize, m_nsize;
  Sizes *m_mlsize, *m_nlsize;
#else
  unsigned int m_numDims;
  SizesArray m_sizes;
  SizesArray m_lsizes;
#endif
  Name m_varName;
 public:
  InputNode();
#if TWOD
  InputNode(NodeType type, Size m, Size n, string name);
#if DODM
  InputNode(NodeType type, Size m, Size n, string name, DistType dist);
#endif
#else
  InputNode(NodeType type, unsigned int numDims, const SizesArray sizes, string name);
  InputNode(NodeType type, unsigned int numDims, const SizesArray sizes, DistType dist, string name);
#endif
  virtual ~InputNode();
  virtual NodeType GetType() const {return m_type;}
  static Node* BlankInst() { return  new InputNode; }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
#if DODM
  virtual DistType GetDistType(unsigned int num) const { return m_varName.m_type; }
#endif
  virtual void SanityCheck();
  virtual void Prop();
  virtual void PrintCode(IndStream &out) {}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "inputNode";}
#if TWOD
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
#else
  virtual const unsigned int NumDims(unsigned int num) const;
  virtual const Sizes* Len(unsigned int num, unsigned int dim) const;
  virtual const Sizes* LocalLen(unsigned int num, unsigned int dim) const;
#endif
  virtual Name GetName(unsigned int num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual void ClearSizeCache();
  virtual void BuildSizeCache();
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
};

class OutputNode : public DLANode
{
  NodeType m_type;
 public:
  OutputNode() : m_type("OutputNode") {}
  OutputNode(NodeType type) : m_type(type) {}
  virtual NodeType GetType() const {return m_type;}
  static Node* BlankInst() { return  new OutputNode; }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
#if DODM
  virtual DistType GetDistType(unsigned int num) const;
#endif
  virtual void SanityCheck();
  virtual void Prop();
  virtual void PrintCode(IndStream &out) {}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "outputNode";}
#if TWOD
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
#else
  virtual const unsigned int NumDims(unsigned int num) const;
  virtual const Sizes* Len(unsigned int num, unsigned int dim) const;
  virtual const Sizes* LocalLen(unsigned int num, unsigned int dim) const;
#endif
  virtual Name GetName(unsigned int num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
};

#if TWOD


//Constant value (e.g. used for Axpy)
class ConstVal : public DLANode
{
  Name m_varName;
  Coef m_val;
 public:
  ConstVal(string name, Coef val);
  virtual NodeType GetType() const {return "const val";}
  static Node* BlankInst() { return  new ConstVal("const", COEFZERO); }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
#if DOELEM
  virtual DistType GetDistType(unsigned int num) const { return D_STAR_STAR; }
#endif
  virtual void SanityCheck();
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "const val";}
#if TWOD
  virtual const Sizes* GetM(unsigned int num) const {return ONES;}
  virtual const Sizes* GetN(unsigned int num) const {return ONES;}
  virtual const Sizes* LocalM(unsigned int num) const {return ONES;}
  virtual const Sizes* LocalN(unsigned int num) const {return ONES;}
#else
blah
#endif
  virtual Name GetName(unsigned int num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
};

class TempVarNode : public DLANode
{
#if DODM
  DistType m_distType;
#endif
  string m_name;
  Sizes *m_mlsize, *m_nlsize;
 public:
 TempVarNode() 
   : 
#if DOELEM
m_distType(D_LASTDIST),
#elif DOTENSORS
  m_distType(DEFAULTDISTTYPE),
#endif
 m_mlsize(NULL), m_nlsize(NULL) {}

 TempVarNode(string name) 
   :  
#if DOELEM
m_distType(D_LASTDIST),
#elif DOTENSORS
  m_distType(DEFAULTDISTTYPE),
#endif
  m_name(name), m_mlsize(NULL), m_nlsize(NULL) {}

#if DODM
 TempVarNode(DistType dist) 
   : m_distType(dist), m_mlsize(NULL), m_nlsize(NULL) {}
 TempVarNode(DistType dist, string name) :  m_distType(dist), m_name(name), m_mlsize(NULL), m_nlsize(NULL) {}
#endif
  virtual NodeType GetType() const;
  static Node* BlankInst() { return  new TempVarNode; }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
#if DODM
  virtual DistType GetDistType(unsigned int num) const { return m_distType; }
#endif
  virtual void SanityCheck();
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "tempVar";}
#if TWOD
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
#else
basfd
#endif
  virtual Name GetName(unsigned int num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual void ClearSizeCache();
  virtual void BuildSizeCache();
  virtual bool IsReadOnly() const {return true;}
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
  virtual bool IsDataDependencyOfInput() const {return false;}
};


#if DOELEM
class MakeTrapNode : public DLAOp<1,1>
{
  Side m_side;
  Tri m_tri;
  int m_offset;
 public:
 MakeTrapNode(Side side, Tri tri, int offset) : m_side(side), m_tri(tri), m_offset(offset) {}
  virtual DistType GetDistType(unsigned int num) const { return InputDistType(0); }
  static Node* BlankInst() { return  new MakeTrapNode(LEFT,LOWER,0); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const {return "MakeTrap";}
  virtual void SanityCheck();
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
  bool CanApply(const Poss *poss, const Node *node) const;
  void Apply(Poss *poss, Node *node) const;
};
#endif

#if DOELEM||DOBLIS
class ScaleNode : public DLAOp<1,1>
{
 public:
  Coef m_val;
  ScaleNode(Layer layer, Coef val);
#if DOELEM
  virtual DistType GetDistType(unsigned int num) const;
#endif
  static Node* BlankInst() { return  new ScaleNode(ABSLAYER, COEFZERO); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual NodeType GetType() const {return "Scale";}
  virtual void SanityCheck();
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
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return false;}
};

class ScaleTrapNode : public DLAOp<1,1>
{
  Side m_side;
  Tri m_tri;
  Coef m_val;
 public:
 ScaleTrapNode(Layer layer, Side side, Tri tri, Coef val) : m_side(side), m_tri(tri), m_val(val) {SetLayer(layer);}
#if DOELEM
  virtual DistType GetDistType(unsigned int num) const;
#endif
  static Node* BlankInst() { return  new ScaleTrapNode(ABSLAYER,LEFT,LOWER,COEFONE); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual NodeType GetType() const {return "ScaleTrap";}
  virtual void SanityCheck();
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
  Sizes *m_sizes,*m_lsizes;
 ViewPan(bool isVert, string name) 
   : m_isVert(isVert), m_name(name), m_sizes(NULL), m_lsizes(NULL) {}
  virtual ~ViewPan();
  virtual void Prop();
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual unsigned int NumOutputs() const {return 1;}
  virtual void ClearSizeCache();
  virtual void BuildSizeCache();
#if TWOD
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
#else
sdlkfj
#endif
  virtual Name GetName(unsigned int num) const;
  virtual NodeType GetType() const {return "ViewPan";}
  static Node* BlankInst() { return new ViewPan(false,""); }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
#if DOELEM
  virtual DistType GetDistType(unsigned int num) const { return InputDistType(0); }
#endif
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ViewPan";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
};


//View the square around the diagonal on block
class ViewAroundDiag : public DLANode
{
 public:
  bool m_isVert;
  string m_name;
  Sizes *m_sizes0,*m_lsizes0;
  Sizes *m_sizes1,*m_lsizes1;
 ViewAroundDiag(bool isVert, string name) 
   : m_isVert(isVert), m_name(name), 
    m_sizes0(NULL), m_lsizes0(NULL),
    m_sizes1(NULL), m_lsizes1(NULL) {}
  virtual ~ViewAroundDiag();
  virtual void Prop();
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual unsigned int NumOutputs() const {return 2;}
  virtual void ClearSizeCache();
  virtual void BuildSizeCache();
#if TWOD
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
#else
sdlkjf
#endif
  virtual Name GetName(unsigned int num) const;
  virtual NodeType GetType() const {return "ViewAround";}
  static Node* BlankInst() { return  new ViewAroundDiag(false,""); }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
#if DOELEM
  virtual DistType GetDistType(unsigned int num) const { return InputDistType(0); }
#endif
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ViewAroundDiag";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
};

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
#if DOELEM
  virtual DistType GetDistType(unsigned int num) const { return InputDistType(2+num); }
#endif
  virtual void PrintCode(IndStream &out) {}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ViewAroundDiagCombine";}
  virtual NodeType GetType() const {return "ViewAroundCombine";}
  virtual void Prop();
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
};
#endif


class ViewTL : public DLANode
{
 public:
  ViewTL(Layer layer) {SetLayer(layer);}
#if TWOD
  virtual const Sizes* GetM(unsigned int num) const
  { return GetInputM(2); }
  virtual const Sizes* GetN(unsigned int num) const
  { return GetInputN(1); }
  virtual const Sizes* LocalM(unsigned int num) const
  { return InputLocalM(2); }
  virtual const Sizes* LocalN(unsigned int num) const
  { return InputLocalN(1); }
#else
bljsdf
#endif
  virtual Name GetName(unsigned int num) const;
  virtual void Prop();
  virtual void SanityCheck();
  virtual unsigned int NumOutputs() const {return 1;}
  static Node* BlankInst() { return  new ViewTL(ABSLAYER); }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
#if DOELEM
  virtual DistType GetDistType(unsigned int num) const { return InputDistType(0); }
#endif
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ViewTL";}
  virtual NodeType GetType() const {return "ViewTL";}
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
};


class ViewTLCombine : public DLANode
{
 public:
  ViewTLCombine(Layer layer) {SetLayer(layer);}

  virtual const Sizes* GetM(unsigned int num) const
  { return GetInputM(1); }
  virtual const Sizes* GetN(unsigned int num) const
  { return GetInputN(1); }
  virtual const Sizes* LocalM(unsigned int num) const
  { return InputLocalM(1); }
  virtual const Sizes* LocalN(unsigned int num) const
  { return InputLocalN(1); }
  virtual Name GetName(unsigned int num) const {return GetInputName(1);}
  virtual void Prop();
  virtual void SanityCheck() {DLANode::SanityCheck();}
  virtual unsigned int NumOutputs() const {return 1;}
  static Node* BlankInst() { return  new ViewTLCombine(ABSLAYER); }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
#if DOELEM
  virtual DistType GetDistType(unsigned int num) const { return InputDistType(1); }
#endif
  virtual void PrintCode(IndStream &out) {}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ViewTLCombine";}
  virtual NodeType GetType() const {return "ViewTLCombine";}
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
};
#endif
