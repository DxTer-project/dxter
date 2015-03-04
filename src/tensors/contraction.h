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
#if DOTENSORS

#include "DLAOp.h"
#include "transform.h"
#include "lowerLayer.h"

class Contraction : public DLAOp<3,1>
{
 public:
  Coef m_alpha, m_beta;
  Type m_type;
  string m_AIndices, m_BIndices, m_CIndices;
  string m_contIndices;
  bool m_needsPacking;
  Contraction(Layer layer, Coef alpha, Coef beta, Type type, 
	      string AIndices, string BIndices, string CIndices,
	      string contIndices);
  static Node* BlankInst();
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Contraction";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual NodeType GetType() const;
  virtual Phase MaxPhase() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  void CheckInputTypesAlign() const;
  virtual void AddVariables(VarSet &set) const;
};

class DistContToLocalContStatASumScatter : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  DistContToLocalContStatASumScatter(Layer fromLayer, Layer toLayer)
    : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
  virtual Cost RHSCostEstimate(const Node *node) const;
};


class DistContToLocalContStatBSumScatter : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  DistContToLocalContStatBSumScatter(Layer fromLayer, Layer toLayer)
    : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
  virtual Cost RHSCostEstimate(const Node *node) const;
};

class DistContToLocalContStatC : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  DistContToLocalContStatC(Layer fromLayer, Layer toLayer)
    : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
  virtual Cost RHSCostEstimate(const Node *node) const;
};

/*
class DistContToLocalContStatC : public VarTrans
{
 public:
  virtual string GetType() const {return "Dist Cont to Local Stat C";}
  virtual int CanApply(const Node *node, void **cache) const;
  virtual void Apply(int num, Node *node, void **cache) const;
  virtual bool IsRef() const {return true;}
  virtual void CleanCache(void **cache) const;
  //  virtual Cost RHSCostEstimate(const Node *node) const;
};
*/


class ContractionLoopExp : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  unsigned int m_dim;
  ContractionLoopExp(Layer fromLayer, Layer toLayer, unsigned int dim);
  
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

class ContractionLowerLayer : public LowerLayer
{
 public:
 ContractionLowerLayer(Layer fromLayer, Layer toLayer)
   : LowerLayer(fromLayer, toLayer, 1) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class PermuteWhileUnpacking : public SingleTrans
{
 public:
  virtual string GetType() const { return "PermuteWhileUnpacking"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};


#endif

