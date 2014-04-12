#include "tensorRedist.h"
#include <ostream>
#include <sstream>
#include <algorithm>
#include "helperNodes.h"

#if DOTENSORS

bool SameDims(DistEntry entry1, DistEntry entry2)
{
  DimVec vec1 = entry1.DistEntryDims();
  DimVec vec2 = entry2.DistEntryDims();
  if (vec1.size() != vec2.size())
    return false;
  DimVecIter iter = vec1.begin();
  for(; iter != vec1.end(); ++iter) {
    bool found = false;
    DimVecIter iter2 = vec2.begin();
    for( ; !found && iter2 != vec2.end(); ++iter2) {
      if (*iter == *iter2)
	found = true;
    }
    if (!found)
      return false;
  }
  return true;
}

/*
  void GetSumScatterInfo(const DistType &inType, const DistType &outType, EntrySet &sumDims)
  {
  for(Dim dim = inType.m_numDims-1; dim >= outType.m_numDims; --dim) {
  sumDims.insert(inType.m_dists[dim]);
  }
  }
*/

SumScatterUpdateNode::SumScatterUpdateNode(Coef coef, const EntrySet &sumDims)
  : DLAOp<2,1>(), m_coef(coef), m_sumDims(sumDims)
{
}


void SumScatterUpdateNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const SumScatterUpdateNode *node = (SumScatterUpdateNode*)orig;
  DLAOp<2,1>::Duplicate(node, shallow, possMerging);
  m_coef = node->m_coef;
  m_sumDims = node->m_sumDims;
}

NodeType SumScatterUpdateNode::GetType() const
{
  return "SumScatterUpdate";
}
void SumScatterUpdateNode::SanityCheck()
{
  throw;
}



void SumScatterUpdateNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    m_cost = 0;

    const DistType &inType = InputDistType(0);
    const DistType &outType = InputDistType(1);
    
    if (m_sumDims.empty())
      throw;

    if (!IsScalar(0)) {
      if (m_sumDims.size() != (inType.m_numDims - outType.m_numDims)) {
	cout << "in: " << inType.str() << endl;
	cout << "out: " << outType.str() << endl;
	cout << m_sumDims.size() << " sum dims\n";
	throw;
      }
    }
    
    if (m_sumDims.size() == 1) {
      unsigned int numProcs = 1;

      DimVec sums = (*(m_sumDims.begin())).DistEntryDims();
      DimVecConstIter iter = sums.begin();
      for(; iter != sums.end(); ++iter) {
	GridLens[*iter];
      }

      DLANode *input = (DLANode*)(Input(1));
      unsigned int num = InputConnNum(1);
      const Dim numDims = input->NumDims(num);

      if (input->IsScalar(num)) {
	input = (DLANode*)(Input(0));
	num = InputConnNum(0);
	m_cost = input->LocalLen(num,0)->NumSizes() * ReduceScatter(numProcs,numProcs);
      }
      else {
	const unsigned int totNumIters = input->LocalLen(num,0)->NumSizes();
	for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
	  Cost temp = 1;
	  for (Dim dim = 0; dim < numDims; ++dim) {
	    temp *= (*(input->LocalLen(num,dim)))[iteration];
	  }
	  m_cost += ReduceScatter(temp*numProcs, numProcs);
	}
      }
    }
    else if (inType.m_numDims == outType.m_numDims) {
      cout << "in : " << inType.str() << endl;
      cout << "out : " << outType.str() << endl;
      cout << "sumDims " << m_sumDims.size() << endl;
      cout << Input(0)->GetNodeClass() << endl;
      SumScatterUpdateNode *in = (SumScatterUpdateNode*)Input(0);
      cout << "sumdims = " << in->m_sumDims.size() << endl;
      cout << in->InputDistType(0).str() << endl;
      cout << in->InputDistType(1).str() << endl;
      throw;
    }
  }
}

Phase SumScatterUpdateNode::MaxPhase() const
{

  const DistType &inType = InputDistType(0);
  const DistType &outType = InputDistType(1);

  if (CurrPhase >= SUMSCATTERTENSORPHASE) {
    DistEntry sumDims;
    if (m_sumDims.size() != 1) {
      /*
      if (!m_hasRefined) {
	cout << "Too many SumScatter dimensions\n";
	cout << "1trying SumScatter from " << DistTypeToStr(inType)
	     << " to " << DistTypeToStr(outType)  << endl;
      }
      */
      return SUMSCATTERTENSORPHASE;
    }
    else if (!IsScalar(0)) {
      sumDims = *(m_sumDims.begin());

      bool foundSumDim = false;
      for (Dim dim = 0; !foundSumDim && dim < inType.m_numDims; ++dim) {
	if (dim < outType.m_numDims) {
	  DistEntry inEntry = inType.m_dists[dim];
	  DistEntry outEntry = outType.m_dists[dim];
	  if (inEntry != outEntry) {
	    if (inEntry.IsStar()) {
	      //[*] -> ...
	      if (outEntry != sumDims) {
		cout << "2trying SumScatter from " << DistTypeToStr(inType)
		     << " to " << DistTypeToStr(outType)  << endl;
		cout << inEntry.DistEntryToStr() << " -> "
		     << outEntry.DistEntryToStr() << endl;
		cout << "sumDims " << sumDims.DistEntryToStr() << endl;
		
		throw;
	      }
	    }
	    else {
	      DimVec inVec = inEntry.DistEntryDims();
	      DimVec dims = sumDims.DistEntryDims();

	      DimVec outVec = outEntry.DistEntryDims();
	      inVec.insert(inVec.end(), dims.begin(), dims.end());
	      if (inVec != outVec) {
		/*
		if (!m_hasRefined) {
		  cout << "3trying SumScatter from " << DistTypeToStr(inType)
		       << " to " << DistTypeToStr(outType)  << endl;
		  cout << m_sumDims.size() << " sum dims\n";
		  if (IsScalar(0)) 
		    cout << "is scalar\n";
		  else
		    cout << "is not scalar\n";
		      
		}
		*/

		return SUMSCATTERTENSORPHASE;	
	      }
	    }
	    foundSumDim = true;
	  }
	}
      }
      if (!foundSumDim)
	throw;
    }
  }

  return NUMPHASES;
}

void SumScatterUpdateNode::PrintCode(IndStream &out)
{
  out.Indent();
  *out << GetInputNameStr(1) << ".SumScatterUpdate( ";
  out << m_coef;
  *out << ", " << GetInputNameStr(0) << ", D";

  EntrySetIter setIter = m_sumDims.begin();
  for(; setIter != m_sumDims.end(); ++setIter) {
    DimVec sumDims = (*setIter).DistEntryDims();
    DimVecConstIter iter = sumDims.begin();
    for(; iter != sumDims.end(); ++iter) {
      *out << "_" << *iter;
    }
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



bool SeparateRedistFromSumScatter::CanApply(const Poss *poss, const Node *node) const
{
  const SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;
  
  const DistType &inType = sum->InputDistType(0);
  const DistType &outType = sum->InputDistType(1);

  //  cout << "eval with in " << inType.str() << " and out " << outType.str() << endl;

  EntrySet sumDims;
  //  GetSumScatterInfo(inType, outType, sumDims);

  if (outType.m_numDims) {

    for(Dim dim = 0; dim < min(inType.m_numDims,outType.m_numDims); ++dim) {
      DistEntry inEntry = inType.m_dists[dim];
      DistEntry outEntry = outType.m_dists[dim];
      if (inEntry != outEntry) {
	if (inEntry.IsStar()) {
	  bool found = false;
	  EntrySetIter iter2 = sumDims.begin();
	  for(; !found && iter2 != sumDims.end(); ++iter2) {
	    //in any order
	    if (SameDims(*iter2, outEntry))
	      found = true;
	  }
	  if (!found)
	    return true;
	}
	else if (outEntry.IsStar()) {
	  return true;
	}
	else {
	  DimVec inVec = inEntry.DistEntryDims();
	  DimSet inSet;
	  inSet.insert(inVec.begin(), inVec.end());
	  DimVec outVec = outEntry.DistEntryDims();
	  DimSet outSet;
	  outSet.insert(outVec.begin(), outVec.end());
	
	  bool found = false;
	  EntrySetIter iter = sumDims.begin();
	  for(; !found && iter != sumDims.end(); ++iter) {
	    DimVec sumVec = (*iter).DistEntryDims();
	    DimSet sumSet;
	    sumSet.insert(sumVec.begin(), sumVec.end());

	    if (includes(outSet.begin(), outSet.end(),
			 sumSet.begin(), sumSet.end()))
	      {
		DimVec tempVec = outVec;
		DimVecIter iter2 = tempVec.begin();
		for( ; iter2 != tempVec.end(); ++iter2)  {
		  if (sumSet.find(*iter2) == sumSet.end()) {
		    if (inSet.find(*iter2) == inSet.end())
		      return true;
		  }
		}
		if (inSet.size() != (outVec.size() + sumSet.size()))
		  return true;
		found = true;
	      }
	  }
	}
      }
    }
  }  
  return false;
}

void SeparateRedistFromSumScatter::Apply(Poss *poss, Node *node) const
{
  SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;
  
  DistType inTypeInt = sum->InputDistType(0);
  const DistType &outType = sum->InputDistType(1);

  const EntrySet &sumDims = sum->m_sumDims;

  for(Dim dim = 0; dim < outType.m_numDims; ++dim) {
    DistEntry inEntry = inTypeInt.m_dists[dim];
    DistEntry outEntry = outType.m_dists[dim];
    if (inEntry != outEntry) {
      if (inEntry.IsStar()) {
	bool found = false;
	EntrySetIter iter2 = sumDims.begin();
	for(; !found && iter2 != sumDims.end(); ++iter2) {
	  if (SameDims(*iter2, outEntry))
	    found = true;
	}
	if (!found) {
	  inTypeInt.m_dists[dim] = outType.m_dists[dim];
	}
      }
      else if (outEntry.IsStar()) {
	inTypeInt.m_dists[dim].SetToStar();
      }
      else {
	DimVec inVec = inEntry.DistEntryDims();
	DimSet inSet;
	inSet.insert(inVec.begin(), inVec.end());
	DimVec outVec = outEntry.DistEntryDims();
	DimSet outSet;
	outSet.insert(outVec.begin(), outVec.end());
	
	bool found = false;
	EntrySetIter iter = sumDims.begin();
	for(; !found && iter != sumDims.end(); ++iter) {
	  DimVec sumVec = (*iter).DistEntryDims();
	  DimSet sumSet;
	  sumSet.insert(sumVec.begin(), sumVec.end());

	  if (includes(outSet.begin(), outSet.end(),
		       sumSet.begin(), sumSet.end()))
	    {
	      DimVec tempVec = outVec;
	      DimVecIter iter2 = tempVec.begin();
	      while (iter2 != tempVec.end()) {
		if (sumSet.find(*iter2) != sumSet.end()) {
		  tempVec.erase(iter2);
		  iter2 = tempVec.begin();
		}
		else
		  ++iter2;
	      }
	      inTypeInt.m_dists[dim].DimsToDistEntry(tempVec);
	      found = true;
	    }
	}
      }
    }
  }
  
  
  RedistNode *redist = new RedistNode(inTypeInt);
  redist->AddInput(sum->Input(0), sum->InputConnNum(0));
  
  SumScatterUpdateNode *sum2 = new SumScatterUpdateNode(sum->m_coef, sumDims);
  sum2->AddInput(redist, 0);
  sum2->AddInput(sum->Input(1), sum->InputConnNum(1));

  poss->AddNode(redist);
  poss->AddNode(sum2);
  sum->RedirectChildren(sum2,0);

  node->m_poss->DeleteChildAndCleanUp(node);
}

bool SplitSumScatter::CanApply(const Poss *poss, const Node *node) const
{
  const SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;
  
  const DistType &inType = sum->InputDistType(0);
  const DistType &outType = sum->InputDistType(1);

  if (sum->IsScalar(0)) {
    if (m_dim >= inType.m_numDims)
      return false;
    else
      return sum->m_sumDims.size() > 1;
  }
  else {
    if (m_dim >= outType.m_numDims)
      return false;

    DistEntry inEntry = inType.m_dists[m_dim];
    DistEntry outEntry = outType.m_dists[m_dim];
  
    if (inEntry == outEntry)
      return false;
    else if (outEntry.IsStar())
      return false;

    const EntrySet &sumDims = sum->m_sumDims;
  


    DimVec inVec = inEntry.DistEntryDims();
    DimSet inSet;
    inSet.insert(inVec.begin(), inVec.end());
    DimVec outVec = outEntry.DistEntryDims();
    DimSet outSet;
    outSet.insert(outVec.begin(), outVec.end());

	
    bool found = false;
    EntrySetIter iter = sumDims.begin();
    for(; !found && iter != sumDims.end(); ++iter) {
      DimVec sumVec = (*iter).DistEntryDims();
      DimSet sumSet;
      sumSet.insert(sumVec.begin(), sumVec.end());

      if (includes(outSet.begin(), outSet.end(),
		   sumSet.begin(), sumSet.end()))
	{
	  if (outSet.size() == (inSet.size() + sumSet.size())) {
	    if (includes(outSet.begin(), outSet.end(),
			 inSet.begin(), inSet.end())) 
	      {

		bool needRedist = false;
	      
		if (inType.m_numDims == outType.m_numDims+1) {
		  DimVec tempVec = inVec;
		  tempVec.insert(tempVec.end(), sumVec.begin(), sumVec.end());
		  DimVecIter outIter = outVec.begin();
		  DimVecIter tempIter = tempVec.begin();
		  for(; !needRedist && outIter != outVec.end(); ++outIter, ++tempIter) {
		    if (*outIter != *tempIter) {
		      //		    cout << "applying1 for " << m_dim << " " << inType.str() << " -> " << outType.str() << endl;
		      return true;
		    }
		  }
		  return false;
		}
		//	      cout << "applying2 for " << m_dim << " " << inType.str() << " -> " << outType.str() << endl;
		return true;
	      }
	  }
	  return false;
	}
    }
    return false;
  }
}

void SplitSumScatter::Apply(Poss *poss, Node *node) const
{
  SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;
  
  const DistType &inType = sum->InputDistType(0);
  const DistType &outType = sum->InputDistType(1);


  if (sum->IsScalar(0)) {
    DistEntry inEntry = inType.m_dists[m_dim];

    EntrySet sums = sum->m_sumDims;
    EntrySetIter sumIter = sums.begin();
    if (!sums.erase(inEntry)) {
      cout << "scalar " << sum->GetNameStr(0) << endl;
      cout << "didn't find entry\n";
      cout << "in " << inType.str()
	   << "\nout " << outType.str()
	   << "\ninEntry: " << inEntry.DistEntryToStr() << endl;
      throw;
    }

    if (sums.empty())
      throw;

    TempVarNode *temp = new TempVarNode(sum->GetDistType(0), sums, sum->GetInputName(0).m_name);
    temp->AddInput(node->Input(1), node->InputConnNum(1));
    poss->AddNode(temp);

    EntrySet set;
    set.insert(inEntry);
    
    SumScatterUpdateNode *newSum = new SumScatterUpdateNode(COEFONE, set);
    newSum->AddInput(node->Input(0), node->InputConnNum(0));
    newSum->AddInput(temp, 0);
    poss->AddNode(newSum);

	 
    node->ChangeInput2Way(node->Input(0), node->InputConnNum(0),
			  newSum, 0);
    sum->m_sumDims = sums;
    return;
  }
  else {
    DistEntry inEntry = inType.m_dists[m_dim];
    DistEntry outEntry = outType.m_dists[m_dim];
  
    DimVec inVec = inEntry.DistEntryDims();
    DimSet inSet;
    inSet.insert(inVec.begin(), inVec.end());
    DimVec outVec = outEntry.DistEntryDims();
    DimSet outSet;
    outSet.insert(outVec.begin(), outVec.end());
	
    EntrySetIter iter = sum->m_sumDims.begin();
    for(; iter != sum->m_sumDims.end(); ++iter) {
      DimVec sumVec = (*iter).DistEntryDims();
      DimSet sumSet;
      sumSet.insert(sumVec.begin(), sumVec.end());

      if (includes(outSet.begin(), outSet.end(),
		   sumSet.begin(), sumSet.end()))
	{
	  DimVec tempVec = inVec;
	  tempVec.insert(tempVec.end(), sumVec.begin(), sumVec.end());
	
	  if (tempVec.size() != outVec.size())
	    throw;
	
	  bool needRedist = false;
	
	  DimVecIter outIter = outVec.begin();
	  DimVecIter tempIter = tempVec.begin();
	  for(; !needRedist && outIter != outVec.end(); ++outIter, ++tempIter) {
	    if (*outIter != *tempIter)
	      needRedist = true;
	  }

	
	  DistType newOutTypeFull;
	  DistType newOutTypeNoSums;
	  newOutTypeFull.PrepForNumDims(inType.m_numDims - 1);
	  newOutTypeNoSums.PrepForNumDims(outType.m_numDims);
	  Dim dimLHS = 0;
	  Dim dimRHS = 0;
	  while (dimRHS < inType.m_numDims) {
	    if (dimRHS == m_dim) {
	      newOutTypeNoSums.m_dists[dimLHS] = outType.m_dists[m_dim];
	      newOutTypeFull.m_dists[dimLHS] = outType.m_dists[m_dim];
	      ++dimLHS;
	      ++dimRHS;
	    }
	    else {
	      if (inType.m_dists[dimRHS] == *iter) {
		++dimRHS;
	      }
	      else {
		if (dimRHS < outType.m_numDims)
		  newOutTypeNoSums.m_dists[dimLHS] = inType.m_dists[dimRHS];
		newOutTypeFull.m_dists[dimLHS] = inType.m_dists[dimRHS];
		++dimLHS;
		++dimRHS;
	      }
	    }
	  }

	  //	cout << "now " << inType.str() << " " << newOutTypeFull.str() << endl;


	  EntrySet sums = sum->m_sumDims;
	  //	cout << "sums before " << sums.size() << endl;
	  //	cout << "sumDims " << sum->m_sumDims.size() << endl;
	  EntrySetIter sumIter = sums.begin();
	  for(; sumIter != sums.end(); ++sumIter) {
	    if (*sumIter == *iter) {
	      sums.erase(sumIter);
	      break;
	    }
	  }
	  //	cout << "sums after " << sums.size() << endl;

	  TempVarNode *temp = new TempVarNode(newOutTypeNoSums, sums, node->GetInputName(0).m_name);
	  temp->AddInput(node->Input(1), node->InputConnNum(1));
	  poss->AddNode(temp);


	  EntrySet set;
	  set.insert(*iter);
	  SumScatterUpdateNode *newSum = new SumScatterUpdateNode(COEFONE, set);
	  newSum->AddInput(node->Input(0), node->InputConnNum(0));
	  newSum->AddInput(temp, 0);
	  poss->AddNode(newSum);

	  RedistNode *newRedist = NULL;
	  if (needRedist) {
	    newRedist = new RedistNode(newOutTypeFull);
	    newRedist->AddInput(newSum, 0);
	    poss->AddNode(newRedist);
	  }

	  if (DistTypeNotEqual(newOutTypeFull, outType)) {
	    if (needRedist)
	      node->ChangeInput2Way(node->Input(0), node->InputConnNum(0),
				    newRedist, 0);
	    else
	      node->ChangeInput2Way(node->Input(0), node->InputConnNum(0),
				    newSum, 0);

	    sum->m_sumDims = sums;
	  }
	  else {
	    if (needRedist)
	      sum->RedirectChildren(newRedist,0);
	    else
	      sum->RedirectChildren(newSum,0);
	  
	    node->m_poss->DeleteChildAndCleanUp(node);
	  }
	  return;
	}
    }
  }
}


#endif





