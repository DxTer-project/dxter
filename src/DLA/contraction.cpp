#include "contraction.h"

#if DOTENSORS

Contraction::Contraction(Layer layer, Coef alpha, Coef beta, Type type, string indices)
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
  return "Contraction " + m_indices
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
    m_cost = 0;
    cout << "improve Contraction::Prop code\n";
    cout << "reflect in DistContToLocalContStatC::RHSCostEstimate\n";
    IndexDimMap dims = MapIndicesToDims(m_indices,GetInputName(0).m_indices);
    const Sizes *sizes = InputLocalLen(2,0);
    unsigned int totNumIters = sizes->NumSizes();
    for(unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
      Cost temp = 1;
      Dim numDims = InputNumDims(2);
      for (Dim dim = 1; dim < numDims; ++dim) {
	temp *= (*InputLocalLen(2,dim))[iteration];
      }
      IndexDimMapConstIter iter = dims.begin();
      for(; iter != dims.end(); ++iter) {
	temp *= (*InputLocalLen(0,*iter))[iteration];
      }
      m_cost += temp;
    }
    m_cost *= 2;
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

bool DistContToLocalContStatC::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != Contraction::GetClass())
    throw;
  const Contraction *cont = (Contraction*)node;
  Dim numContDims = cont->m_indices.length();
  if (numContDims > (MAX_NUM_DIMS / 4))
    throw;

  
  
  sdlfkj
}

void DistContToLocalContStatC::Apply(Poss *poss, Node *node) const
{
  adsflkj
}

Cost DistContToLocalContStatC::RHSCostEstimate(const Node *node) const
{
  Cost cost = 0;
  const Contraction *cont = (Contraction*)node;
  IndexDimMap dims = MapIndicesToDims(cont->m_indices,cont->GetInputName(0).m_indices);
  const Sizes *sizes = cont->InputLocalLen(2,0);
  unsigned int totNumIters = sizes->NumSizes();
  for(unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
    Cost temp = 1;
    Dim numDims = cont->InputNumDims(2);
    for (Dim dim = 1; dim < numDims; ++dim) {
      temp *= (*(cont->InputLocalLen(2,dim)))[iteration];
    }
    IndexDimMapConstIter iter = dims.begin();
    for(; iter != dims.end(); ++iter) {
      temp *= (*(cont->InputLocalLen(0,*iter)))[iteration];
    }
    cost += temp;
  }
  return cost * 2;
}


#endif
