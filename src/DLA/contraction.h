#pragma once

#include "layers.h"
#if DOTENSORS

#include "DLAOp.h"
#include "transform.h"

class Contraction : public DLAOp<3,1>
{
 public:
  Coef m_alpha, m_beta;
  Type m_type;
  string m_indices;
  Contraction(Layer layer, Coef alpha, Coef beta, Type type, string indices);
  static Node* BlankInst();
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Contraction";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual NodeType GetType() const;
  virtual const DistType& GetDistType(unsigned int num) const;
  virtual Phase MaxPhase() const;
  virtual void SanityCheck();
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
};

class DistContToLocalContStatAAllReduce : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  DistContToLocalContStatAAllReduce(Layer fromLayer, Layer toLayer)
    : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

class DistContToLocalContStatASumScatter : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  DistContToLocalContStatASumScatter(Layer fromLayer, Layer toLayer)
    : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

class DistContToLocalContStatC : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  DistContToLocalContStatC(Layer fromLayer, Layer toLayer)
    : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

/*
class DistContToLocalContStatC : public VarTrans
{
 public:
  virtual string GetType() const {return "Dist Cont to Local Stat C";}
  virtual int CanApply(const Poss *poss, const Node *node, void **cache) const;
  virtual void Apply(Poss *poss, int num, Node *node, void **cache) const;
  virtual bool IsRef() const {return true;}
  virtual void CleanCache(void **cache) const;
  //  virtual Cost RHSCostEstimate(const Node *node) const;
};
*/
#endif

