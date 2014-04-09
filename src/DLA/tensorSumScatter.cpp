#include "tensorRedist.h"
#include <ostream>
#include <sstream>

#if DOTENSORS

SumScatterUpdateNode::SumScatterUpdateNode(const DimSet &sumDims, Coef coef)
  : DLAOp<2,1>(), m_coef(coef)
{
  m_sumDims = sumDims;
}


void SumScatterUpdateNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const SumScatterUpdateNode *node = (SumScatterUpdateNode*)orig;
  DLAOp<2,1>::Duplicate(node, shallow, possMerging);
  m_sumDims = node->m_sumDims;  
  m_coef = node->m_coef;
}

NodeType SumScatterUpdateNode::GetType() const
{
  stringstream str;
  str << "SumScatterUpdate";
  DimSetConstIter iter = m_sumDims.begin();
  for(; iter != m_sumDims.end(); ++iter)
    str << *iter << ",";
  return str.str();
}
void SumScatterUpdateNode::SanityCheck()
{
  throw;
}

void SumScatterUpdateNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    const DistType &inType = InputDistType(0);
    const DistType &outType = InputDistType(1);

    if (inType.m_numDims != (outType.m_numDims + m_sumDims.size())) {
      cout << "inType " << inType.str() << endl;
      cout << "outType " << outType.str() << endl;
      cout << "sumDims ";
      DimSetIter iter = m_sumDims.begin();
      for( ; iter != m_sumDims.end(); ++iter) 
	cout << *iter << " ";
      cout << endl;      
      throw;
    }

#if 1
    if (CurrPhase >= ROTENSORPHASE) {
      for (Dim dim = 0; dim < inType.m_numDims; ++dim) {
	if (inType.m_dists[dim] != outType.m_dists[dim]) {
	  if (inType.m_dists[dim] != 0) {
	    cout << "trying SumScatter from " << DistTypeToStr(inType)
		 << " to " << DistTypeToStr(outType) << " with reduction on ";
	    DimSetIter iter = m_sumDims.begin();
	    for( ; iter != m_sumDims.end(); ++iter)
	      cout << *iter << " ";
	    cout << endl;
	    throw;
	  }
	  else {
	    DimVec dims = DistType::DistEntryDims(outType.m_dists[dim]);
	    DimVecIter iter = dims.begin();
	    for(; iter != dims.end(); ++iter) {
	      if (m_sumDims.find(*iter) == m_sumDims.end()) {
		cout << "SumScatter error\n";
		cout << "trying SumScatter from " << DistTypeToStr(inType)
		     << " to " << DistTypeToStr(outType) << " with reduction on ";
		DimSetIter iter = m_sumDims.begin();
		for( ; iter != m_sumDims.end(); ++iter)
		  cout << *iter << " ";
		cout << endl;
		throw;
	      }
	    }
	  }
	}
      }
    }
#endif


    m_cost = 0;
    unsigned int numProcs = 1;

    DimSetIter iter = m_sumDims.begin();
    for(; iter != m_sumDims.end(); ++iter) {
      numProcs *= GridLens[*iter];
    }
      
    DLANode *input = (DLANode*)(Input(1));
    unsigned int num = InputConnNum(1);

    const unsigned int totNumIters = input->LocalLen(num,0)->NumSizes();
    const Dim numDims = input->NumDims(num);

    for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
      Cost temp = 1;
      for (Dim dim = 0; dim < numDims; ++dim) {
	temp *= (*(input->LocalLen(num,dim)))[iteration];
      }
      m_cost += AllReduce(temp*numProcs, numProcs);
    }
  }
}

void SumScatterUpdateNode::PrintCode(IndStream &out)
{
  out.Indent();
  *out << GetInputNameStr(1) << ".SumScatterUpdate( ";
  out << m_coef;
  *out << ", " << GetInputNameStr(0)
       << ", m";
  DimSetConstIter iter = m_sumDims.begin();
  for(; iter != m_sumDims.end(); ++iter) {
    *out << "_" << *iter;
  }
  *out << " );\n";
}

void SumScatterUpdateNode::FlattenCore(ofstream &out) const
{
  throw;
}

void SumScatterUpdateNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  throw;
}



#endif
