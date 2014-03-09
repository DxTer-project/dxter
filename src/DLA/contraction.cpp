#include "contraction.h"

#if DOTENSORS

Contraction::Contraction(Layer layer, Coef alpha, Coef beta, Type type, string indiced)
:
  m_alpha(alpha),
  m_beta (beta),
  m_type (type),
  m_indices(indices)
{
  SetLayer(layer);
}

Node* Contraction::BlankInst()
{
  return new Contraction(ABSLAYER, COEFONE, COEFONE, REAL, "");
}

NodeType Contraction::GetType() const
{
  return "Contraction " + m_indices +
    + LayerNumToStr(GetLayer());
}

void Contraction::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const Contraction *cont = (Contraction*)orig;
  m_alpha = cont->m_alpha;
  m_beta = cont->m_beta;
  m_type = cont->m_type;
  m_indices = cont->m_indices;
}

void Contraction::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_alpha);
  WRITE(m_beta);
  WRITE(m_type);
  out << m_indices << endl;
}

void Contraction::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  READ(m_alpha);
  READ(m_beta);
  READ(m_type);
  getline(in, m_indices);
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
    Sizes sizes = InputLocalLen(2,0);
    unsigned int numDims = InputNumDims(2);
    for (unsigned int dim = 1; dim < numDims; ++dim) {
      sizes.MultBy(InputLocalLen(2,dim));
    }
    vector<unsigned int> dims = MapIndicesToDims(m_indices,GetInputName(0).m_indices);
    vector<unsigned int>::const_iterator iter = dims.begin();
    for(; iter != dims.end(); ++iter) {
      sizes.MultBy(InputLocalLen(0,*iter));
    }
    m_cost = 2 * sizes.Sum();
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
