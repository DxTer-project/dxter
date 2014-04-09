#include "tensorRedist.h"
#include <ostream>
#include <sstream>

#if DOTENSORS

SumScatterUpdateNode::SumScatterUpdateNode(Coef coef)
  : DLAOp<2,1>(), m_coef(coef)
{
}


void SumScatterUpdateNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const SumScatterUpdateNode *node = (SumScatterUpdateNode*)orig;
  DLAOp<2,1>::Duplicate(node, shallow, possMerging);
  m_coef = node->m_coef;
}

NodeType SumScatterUpdateNode::GetType() const
{
  return "SumScatterUpdate";
}
void SumScatterUpdateNode::SanityCheck()
{
  throw;
}

Dim GetSumDim(const DistType &inType, const DistType &outType)
{
  if (inType.m_numDims != (outType.m_numDims+1))
    throw;
  for (Dim dim = outType.m_numDims-1; dim > 0; -- dim) {
    unsigned int inEntry = inType.m_dists[dim+1];
    unsigned int outEntry = outType.m_dists[dim];
    if (inEntry != outEntry) {
      if (!inEntry || inEntry > NUM_GRID_DIMS) {
	cout << "trying SumScatter from " << DistTypeToStr(inType)
	     << " to " << DistTypeToStr(outType)  << endl;
	throw;
      }
      else {
	if (inType.m_dists[dim] != outEntry){
	  cout << "trying SumScatter from " << DistTypeToStr(inType)
	       << " to " << DistTypeToStr(outType)  << endl;
	  throw;
	}
	return inEntry - 1;
      }
    }
  }
  cout << "trying SumScatter from " << DistTypeToStr(inType)
       << " to " << DistTypeToStr(outType)  << endl;
  throw;
}

void SumScatterUpdateNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    m_cost = 0;

    const DistType &inType = InputDistType(0);
    const DistType &outType = InputDistType(1);
    
    if (inType.m_numDims == (outType.m_numDims+1)) {
      unsigned int numProcs = GridLens[GetSumDim(inType, outType)];
      
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
}

Phase SumScatterUpdateNode::MaxPhase() const
{

    const DistType &inType = InputDistType(0);
    const DistType &outType = InputDistType(1);

    if (CurrPhase >= ROTENSORPHASE) {
      Dim sumDim;
      if (inType.m_numDims != (outType.m_numDims + 1)) {
	if (!m_hasRefined) {
	  cout << "Too many SumScatter dimensions\n";
	  cout << "trying SumScatter from " << DistTypeToStr(inType)
	       << " to " << DistTypeToStr(outType)  << endl;
	}
	return ROTENSORPHASE;
      }
      else {
	sumDim = GetSumDim(inType, outType);
      }
      bool foundSumDim = false;
      for (Dim dim = 0; dim < inType.m_numDims; ++dim) {
	if (dim < outType.m_numDims) {
	  unsigned int inEntry = inType.m_dists[dim];
	  unsigned int outEntry = outType.m_dists[dim];
	  if (inEntry != outEntry) {
	    if (foundSumDim) {
	      cout << "trying SumScatter from " << DistTypeToStr(inType)
		   << " to " << DistTypeToStr(outType)  << endl;

	      throw;
	    }
	    if (!inEntry) {
	      //[*] -> ...
	      if (outEntry != sumDim) {
		cout << "trying SumScatter from " << DistTypeToStr(inType)
		     << " to " << DistTypeToStr(outType)  << endl;
		
		throw;
	      }
	    }
	    else {
	      DimVec inVec = DistType::DistEntryDims(inEntry);
	      DimVec outVec = DistType::DistEntryDims(outEntry);
	      inVec.push_back(sumDim);
	      if (inVec != outVec) {
		cout << "trying SumScatter from " << DistTypeToStr(inType)
		     << " to " << DistTypeToStr(outType)  << endl;

		throw;
	      }
	    }
	    foundSumDim = true;
	  }
	}
	else {
	  if (!inType.m_dists[dim] || inType.m_dists[dim] > NUM_GRID_DIMS) {
	    throw;
	  }
	}
	
      }
      if (!foundSumDim)
	throw;
    }



  return NUMPHASES;
}

void SumScatterUpdateNode::PrintCode(IndStream &out)
{
  out.Indent();
  *out << GetInputNameStr(1) << ".SumScatterUpdate( ";
  out << m_coef;
  *out << ", " << GetInputNameStr(0);

  const DistType &inType = InputDistType(0);
  const DistType &outType = InputDistType(1);
  if (inType.m_numDims == (outType.m_numDims + 1)) {
    *out << ", D_" << GetSumDim(inType, outType);
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

