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

#include "transform.h"
#include "DLAOp.h"
#include "distributions.h"
#include "gemm.h"
#include "herk.h"
#include "hemm.h"
#include "TrProps.h"

bool IsDMTrxm(const Node *node);
bool IsDMTrxm(const Node *node, bool invert);

class Trxm : public DLAOp<2,1>, public TrProps
{
 public:
  Trxm(bool invert, Layer layer, Side side, Tri tri, Diag diag, Trans trans, Coef coeff, Type type);
  static Node* BlankInst() { return new Trxm(false, ABSLAYER, LEFT,LOWER,UNIT,NORMAL,COEFONE,REAL); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Trxm";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual DistType GetDistType(unsigned int num) const;
  virtual Phase MaxPhase() const;
  virtual bool DoNotCullDP() const;
  virtual NodeType GetType() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual void SanityCheck();

  static Cost GetCost(Layer layer, Side side, const Sizes *localMs, const Sizes *localNs);
  virtual bool ShouldCullSR() const;
};

class Trmm3 : public DLAOp<3,1>, public TrProps
{
 public:
  Coef m_beta;
  Trmm3(Layer layer, Side side, Tri tri, Diag diag, Trans trans, Coef coeff, Coef beta, Type type);
  static Node* BlankInst() { return new Trmm3(ABSLAYER, LEFT,LOWER,UNIT,NORMAL,COEFONE,COEFONE,REAL); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Trmm3";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual DistType GetDistType(unsigned int num) const;
  virtual Phase MaxPhase() const;
  virtual bool DoNotCullDP() const;
  virtual NodeType GetType() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual bool ShouldCullSR() const;
};

class TrxmBP : public DLAOp<3,2>, public TrProps
{
 public:
  Coef m_beta;
  TrxmBP(bool invert, Layer layer, Side side, Tri tri, Trans trans, 
	 Coef coeff, Coef beta, Type type);
  static Node* BlankInst() { return new TrxmBP(false,SQ2LAYER,LEFT,LOWER,NORMAL,COEFONE,COEFONE,REAL); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "TrxmBP";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual DistType GetDistType(unsigned int num) const;
  virtual Phase MaxPhase() const {return NUMPHASES;}
  virtual NodeType GetType() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
};

class TrxmLoopExp : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  unsigned int m_dim;
 TrxmLoopExp(Layer fromLayer, Layer toLayer, unsigned int dim) 
   : m_fromLayer(fromLayer), m_toLayer(toLayer), m_dim(dim) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

class Trmm3LoopExp : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  unsigned int m_dim;
 Trmm3LoopExp(Layer fromLayer, Layer toLayer, unsigned int dim) 
   : m_fromLayer(fromLayer), m_toLayer(toLayer), m_dim(dim) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

class DistTrxmToLocalTrxm : public SingleTrans
{
  DistType m_leftType, m_rightType;
 public:
  DistTrxmToLocalTrxm(DistType leftType, DistType rightType) : m_leftType(leftType), m_rightType(rightType) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
  virtual Cost RHSCostEstimate(const Node *node) const;
};


//Uses special Trsm not in Elemental but something Jack built on top of Elemental
// for large triangular matrix
class DistTrsmToSpecialLocalTrsm : public SingleTrans
{
 public:
  virtual string GetType() const  {return "dist trsm to special local"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};


class LocalTrmmAcc : public DLAOp<3,1>
{
 public:
  Side m_side;
  Tri m_tri;
  Diag m_diag;
  Trans m_trans;
  Coef m_coeff;
  Type m_type;
  LocalTrmmAcc(Side side, Tri tri, Diag diag, Trans trans, Coef coeff, Type type);
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const { return "LocalTrmmAcc"; }
  static Node* BlankInst() { return  new LocalTrmmAcc(LEFT, UPPER, UNIT, TRANS, COEFONE, REAL);} 
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual DistType GetDistType(unsigned int num) const { return InputDistType(2+num); }
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "localTrmmAcc";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  static Cost GetCost(Side side, bool isTrans, const Sizes *localMs, const Sizes *localNs);
};

class DistTrmmToLocalTrmmStatA : public  SingleTrans
{
 public:
  virtual string GetType() const {return "Distributed Trmm to Local Trmm stat A";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
  virtual Cost RHSCostEstimate(const Node *node) const;
};

class TrxmTrans : public TransTransformation
{
 public:
  TrxmTrans(Trans trans) : TransTransformation(1,trans) {}
  virtual string GetTransType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void PreApply(Node *node) const;
  virtual void PostApply(Node *node) const;
};

//Uses Trsm instead of Trmm with triangular matrix inversion
class DTrmmToTrsm : public SingleTrans
{
 public:
  virtual string GetType() const {return "Dist Trmm of inverse to Trsm";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;  
};

class LTrmmToTrsm : public SingleTrans
{
 public:
  virtual string GetType() const {return "Local Trmm of inverse to Trsm";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;  
};

class BLISTrxmLoopExp : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  BLISTrxmLoopExp(Layer fromLayer, Layer toLayer)
    : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const { return "BLISTrxmLoop Exp"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

class BLISTrmm3LoopExp : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  BLISTrmm3LoopExp(Layer fromLayer, Layer toLayer)
    : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const { return "BLISTrmm3Loop Exp"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

template<class TrxmType>
class TrxmRightToLeft : public SingleTrans
{
 public:
  Layer m_layer;
  TrxmRightToLeft(Layer layer)
    : m_layer(layer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class Trmm3RightToLeft : public SingleTrans
{
 public:
  Layer m_layer;
  Trmm3RightToLeft(Layer layer)
    : m_layer(layer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

template<class TrxmType>
class TrxmLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
 TrxmLowerLayer(Layer fromLayer, Layer toLayer)
   : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const { return "Trxm lower layer"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};


class TrmmAxpytoTrxm3 : public SingleTrans
{
 public:
  Layer m_layer;
  TrmmAxpytoTrxm3(Layer layer)
    : m_layer(layer) {}
  virtual string GetType() const { return "TrmmAxpytoTrxm3"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};


class CopyTrmmtoTrxm3 : public SingleTrans
{
 public:
  Layer m_layer;
  CopyTrmmtoTrxm3(Layer layer)
    : m_layer(layer) {}
  virtual string GetType() const { return "CopyTrmmtoTrxm3"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};
