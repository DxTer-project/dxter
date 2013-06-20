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
#include "MPI.h"
#include "gemm.h"
#include "trxm.h"
#include "herk.h"
#include "her2k.h"
#include "hemm.h"
#include "hetrmm.h"

class Axpy : public DLAOp<2,1>
{
 public:
  Coef m_coeff;
 Axpy(Layer layer, Coef coeff) : m_coeff(coeff) { SetLayer(layer); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual NodeType GetType() const;
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Axpy";}
  static Node* BlankInst() { return new Axpy(ABSLAYER, COEFZERO); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual DistType GetDistType(unsigned int num) const;
  virtual void SanityCheck();
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual Phase MaxPhase() const;
  virtual bool ShouldCullDP() const;
  virtual bool DoNotCullDP() const;
  static Cost GetCost(Layer layer, const Sizes *localMs, const Sizes *localNs);
};

class AxpyLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
 AxpyLowerLayer(Layer fromLayer, Layer toLayer)
   : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const { return "Axpy lower layer"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};


class DistAxpyToLocalAxpy : public SingleTrans
{
  DistType m_type;
 public:
  DistAxpyToLocalAxpy(DistType type) : m_type(type) {}
  virtual string GetType() const {return "DistAxpy to LocalAxpy " + DistTypeToStr(m_type);}
  virtual bool WorthApplying(const Node *node) const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
  Cost RHSCostEstimate(const Node *node) const;
};

class AxpyToBLASAxpy : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
 AxpyToBLASAxpy(Layer fromLayer, Layer toLayer) 
   : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const {return "AxpyToBLASAxpy"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};


class Scal : public DLAOp<2,1>
{
 public:
  Scal(Layer layer) {SetLayer(layer);}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Scal";}
  virtual void SanityCheck();
  virtual NodeType GetType() const;
  static Node* BlankInst() { return  new Scal(ABSLAYER); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual DistType GetDistType(unsigned int num) const;
  virtual Phase MaxPhase() const;
  virtual bool ShouldCullDP() const;
};

#if DODPPHASE
class DistScalToLocalScal : public SingleTrans
{
  DistType m_type;
 public:
 DistScalToLocalScal(DistType type) : m_type(type) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};
#endif

class ConstScal : public DLAOp<1,1>
{
 public:
  Coef m_alpha;
 ConstScal(Layer layer, Coef alpha) : m_alpha(alpha) {SetLayer(layer);}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "ConstScal";}
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool ShouldCullDP() const;
  virtual Phase MaxPhase() const;
  virtual NodeType GetType() const;
  static Node* BlankInst() { return new ConstScal(ABSLAYER, COEFZERO); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual DistType GetDistType(unsigned int num) const;
  virtual void SanityCheck();
  virtual void PrintCode(IndStream &out);
  virtual void Prop();
};

#if DODPPHASE
class DistConstScalToLocalConstScal : public SingleTrans
{
  DistType m_type;
 public:
 DistConstScalToLocalConstScal(DistType type) : m_type(type) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};
#endif