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
  static Cost GetCost(Layer layer, const Sizes *localDim1, const Sizes *localDim2, const Sizes *localDim3);
};
#endif
