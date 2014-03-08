#include "contraction.h"

#if DOTENSORS

Contraction::Contraction(Layer layer, Coef alpha, Coef beta, Type type)
:
  m_alpha(alpha),
  m_beta (beta),
  m_type (type)
{
  SetLayer(layer);
}

Node* Contraction::BlankInst()
{
  return new Contraction(ABSLAYER, COEFONE, COEFONE, REAL);
}

NodeType Contraction::GetType() const
{
  return "Contraction "
    + LayerNumToStr(GetLayer());
}

void Contraction::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const Contraction *cont = (Contraction*)orig;
  m_alpha = cont->m_alpha;
  m_beta = cont->m_beta;
  m_type = cont->m_type;
}

void Contraction::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_alpha);
  WRITE(m_beta);
  WRITE(m_type);
}

void Contraction::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  READ(m_alpha);
  READ(m_beta);
  READ(m_type);
}

const DistType& Contraction::GetDistType(unsigned int num) const
{
  return InputDistType(2);
}

void Contraction::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();
    throw;
    //calc cost
  }
}

void Contraction::SanityCheck()
{
  throw;
}

void Contraction::PrintCode(IndStream &out)
{
  throw;
}

Phase Contraction::MaxPhase() const
{
  throw;
}


#endif
