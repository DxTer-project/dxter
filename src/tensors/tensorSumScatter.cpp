#include "tensorRedist.h"
#include <ostream>
#include <sstream>
#include <algorithm>
#include "helperNodes.h"
#include "possTunnel.h"

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

SumScatterUpdateNode::SumScatterUpdateNode(Coef coef, const EntrySet &sumDims)
  : DLAOp<2,1>(), m_coef(coef), m_sumDims(sumDims)
{
  if (!m_sumDims.size())
    throw;
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

    if (inType.m_notReped == outType.m_notReped) {
      if (m_sumDims.size() != (inType.m_numDims - outType.m_numDims)) {
	cout << "in: " << inType.str() << endl;
	cout << "out: " << outType.str() << endl;
	cout << m_sumDims.size() << " sum dims\n";
	throw;
      }



      DimSet allSumDims;
      EntrySetConstIter iter = m_sumDims.begin();
      for(; iter != m_sumDims.end(); ++iter) {
	DimVec vec = (*iter).DistEntryDims();
	allSumDims.insert(vec.begin(), vec.end());
      }
      
      bool foundAll = false;

      for(Dim dim = 0; !foundAll && dim < outType.m_numDims; ++dim) {
	DistEntry outEntry = outType.m_dists[dim];
	DimVec outVec = outEntry.DistEntryDims();
	DimVecIter iter = outVec.begin();
	while (!foundAll && iter != outVec.end()) {
	  DimSetIter found = allSumDims.find(*iter);
	  if (found != allSumDims.end()) {
	    allSumDims.erase(found);
	    iter = outVec.begin();
	    foundAll = allSumDims.empty();
	  }
	  else
	    ++iter;
	}
      }

      if (!foundAll) {
	cout << "Bad sum scatter from " << inType.str() << " -> "
	     << outType.str() << " with " << m_sumDims.size() << " sums\n";
	m_poss->PrintTransVec();
	throw;
      }
    }
    else {
      DistEntry inNotReped = inType.m_notReped;
      DimSet notReped = inNotReped.DistEntryDimSet();
      EntrySetIter iter = m_sumDims.begin();
      for(; iter != m_sumDims.end(); ++iter) {
	DistEntry entry = *iter;
	DimSet set = entry.DistEntryDimSet();
	notReped.insert(set.begin(), set.end());
      }
      DimVec vec;
      vec.insert(vec.end(), notReped.begin(), notReped.end());
      inNotReped.DimsToDistEntry(vec);
      if (outType.m_notReped != inNotReped) {
	cout << "bad notReped\n";
	cout << inType.PrettyStr() << " -> "
	     << outType.PrettyStr() << " with "
	     << m_sumDims.size() << endl;	  
	throw;
      }
      
    }
    
    if (m_sumDims.size() == 1) {
      unsigned int numProcs = 1;

      DimVec sums = (*(m_sumDims.begin())).DistEntryDims();
      DimVecConstIter iter = sums.begin();
      for(; iter != sums.end(); ++iter) {
	numProcs *= GridLens[*iter];
      }

      DLANode *input = (DLANode*)(Input(1));
      unsigned int num = InputConnNum(1);
      const Dim numDims = input->NumDims(num);

      if (inType.m_notReped != outType.m_notReped) {
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

void SumScatterUpdateNode::CheckSumDimsInOutput() const
{
    const DistType &outType = InputDistType(1);
    
    if (InputDistType(0).m_notReped != outType.m_notReped)
      return;

    DimSet allSumDims;
    EntrySetConstIter iter = m_sumDims.begin();
    for(; iter != m_sumDims.end(); ++iter) {
      DimVec vec = (*iter).DistEntryDims();
      allSumDims.insert(vec.begin(), vec.end());
    }
      
    bool foundAll = false;

    for(Dim dim = 0; !foundAll && dim < outType.m_numDims; ++dim) {
      DistEntry outEntry = outType.m_dists[dim];
      DimVec outVec = outEntry.DistEntryDims();
      DimVecIter iter = outVec.begin();
      while (!foundAll && iter != outVec.end()) {
	DimSetIter found = allSumDims.find(*iter);
	if (found != allSumDims.end()) {
	  allSumDims.erase(found);
	  iter = outVec.begin();
	  foundAll = allSumDims.empty();
	}
	else
	  ++iter;
      }
    }
    if (!foundAll) {
      cout << "bad out type\n";
      cout << "outType " << outType.str() << endl;
      cout << "inType " << InputDistType(0).str() << endl;
      m_poss->PrintTransVec();
      throw;
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
    else if (inType.m_notReped == outType.m_notReped) {
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
		cout << inEntry.str() << " -> "
		     << outEntry.str() << endl;
		cout << "sumDims " << sumDims.str() << endl;
		
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
  if (m_sumDims.size() != 1)
    throw;
  out.Indent();

  const DistType &m_srcType = InputDistType(0);
  const DistType &m_destType = InputDistType(1);

  Dim srcNumDims = m_srcType.m_numDims;

  string inName = GetInputName(0).str();
  string outName = GetName(0).str();
  DistEntry sumDims = *(m_sumDims.begin());
    

  *out << "   // " << GetName(0).PrettyStr() 
       << " <- " << GetInputName(0).PrettyStr() 
       << " (with SumScatter on " << sumDims.PrettyStr() << ")\n";

  out.Indent();

  *out << outName << ".SizeTo( " << inName << " );\n";
  out.Indent();
  
  if (m_srcType.m_notReped != m_destType.m_notReped) {
    *out << "ReduceScatterRedistDisappear( " << outName << ", " << inName << ", ";
    bool found = false;
    for(Dim dim = 0; !found &&dim < srcNumDims; ++dim) {
      if (m_srcType.m_dists[dim] == sumDims) {
	*out << dim;
	found = true;
      }
    }
    if (!found)
      throw;
    *out << " );\n";
  }
  else {
    if (srcNumDims != (m_destType.m_numDims+1))
      throw;

    *out << "ReduceScatterRedist( " << outName << ", " << inName << ", ";

    DimVec scatDims = sumDims.DistEntryDims();


    Dim redDim = srcNumDims;
    for (Dim dim = srcNumDims-1; dim >= 0; --dim) {
      if (m_srcType.m_dists[dim] == sumDims) {
	redDim = dim;
	break;
      }
    }
    if (redDim == srcNumDims)
      throw;

    Dim scatDim = srcNumDims;
    for(Dim dim = 0; dim < m_destType.m_numDims; ++dim) {
      DimVec destDims = m_destType.m_dists[dim].DistEntryDims();
      DimVec srcDims = m_srcType.m_dists[dim].DistEntryDims();
      srcDims.insert(srcDims.end(), scatDims.begin(), scatDims.end());
      if (destDims == srcDims) {
	scatDim = dim;
	break;
      }
    }
    if (scatDim == srcNumDims)
      throw;

    *out << redDim << ", " << scatDim << " );\n";
  }
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

  if (inType.m_notReped != outType.m_notReped)
    return false;
  

  const EntrySet &sumDims = sum->m_sumDims;

  DimSet allSumDims;
  EntrySetConstIter iter = sumDims.begin();
  for(; iter != sumDims.end(); ++iter) {
    DimVec vec = (*iter).DistEntryDims();
    allSumDims.insert(vec.begin(), vec.end());
  }

  for(Dim dim = 0; dim < inType.m_numDims; ++dim) {
    DistEntry inEntry = inType.m_dists[dim];
    if (dim >= outType.m_numDims) {
      cout << "dim >= outType.m_numDims\n";
      cout << sum->InputDistType(0).str() << " -> " 
	   << outType.str() << " with " << sum->m_sumDims.size() << endl;
      throw;
    }
    DistEntry outEntry = outType.m_dists[dim];
    if (inEntry != outEntry) {
      DimVec outVec = outEntry.DistEntryDims();
      DimVecIter iter = outVec.begin();
      while (iter != outVec.end()) {
	DimSetIter found = allSumDims.find(*iter);
	if (found != allSumDims.end()) {
	  outVec.erase(iter);
	  allSumDims.erase(found);
	  iter = outVec.begin();
	}
	else
	  ++iter;
      }
      DistEntry entry;
      entry.DimsToDistEntry(outVec);
      if (inType.m_dists[dim] != entry)
	return true;
      if (allSumDims.empty())
	break;
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

  DimSet allSumDims;
  EntrySetConstIter iter = sumDims.begin();
  for(; iter != sumDims.end(); ++iter) {
    DimVec vec = (*iter).DistEntryDims();
    allSumDims.insert(vec.begin(), vec.end());
  }

  for(Dim dim = 0; dim < inTypeInt.m_numDims; ++dim) {
    DistEntry inEntry = inTypeInt.m_dists[dim];
    if (dim >= outType.m_numDims) {
      cout << "dim >= outType.m_numDims\n";
      cout << sum->InputDistType(0).str() << " -> " 
	   << outType.str() << " with " << sum->m_sumDims.size() << endl;
      throw;
    }
    DistEntry outEntry = outType.m_dists[dim];
    if (inEntry != outEntry) {
      DimVec outVec = outEntry.DistEntryDims();
      DimVecIter iter = outVec.begin();
      while (iter != outVec.end()) {
	DimSetIter found = allSumDims.find(*iter);
	if (found != allSumDims.end()) {
	  outVec.erase(iter);
	  allSumDims.erase(found);
	  iter = outVec.begin();
	}
	else
	  ++iter;
      }
      inTypeInt.m_dists[dim].DimsToDistEntry(outVec);
      if (allSumDims.empty())
	break;
    }
  }

  if (sum->Input(1)->m_poss != poss) {
    cout << "input 1 on " << sum->Input(1)->m_poss << " and this is on " << poss << endl;
    throw;
  }
  if (sum->m_poss != poss) {
    cout << "one\n";
    throw;
  }
  
  if (DistTypeEqual(inTypeInt,sum->InputDistType(0))) {
    cout << "same dist\n";
    throw;
  }


  RedistNode *redist = new RedistNode(inTypeInt);
  redist->AddInput(sum->Input(0), sum->InputConnNum(0));
  
  poss->AddNode(redist);

  if (!inTypeInt.IsSane()) {
    cout << "!sane: " << inTypeInt.str() << endl;
    cout << "was " << sum->InputDistType(0).str() << " -> "
	 << outType.str() << " with " << sumDims.size() << endl;
      EntrySetConstIter iter = sumDims.begin();
      for(; iter != sumDims.end(); ++iter) {
	DimVec vec = (*iter).DistEntryDims();
	cout << "sum on " << *(vec.begin()) << endl;
      }

    throw;
  }
  
  SumScatterUpdateNode *sum2 = new SumScatterUpdateNode(sum->m_coef, sumDims);
  sum2->AddInput(redist, 0);
  sum2->AddInput(sum->Input(1), sum->InputConnNum(1));
  poss->AddNode(sum2);


  sum->RedirectChildren(sum2,0);

  node->m_poss->DeleteChildAndCleanUp(node);
}

bool SplitSumScatter::CanApply(const Poss *poss, const Node *node) const
{
  const SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;
  
  const DistType &inType = sum->InputDistType(0);
  const DistType &outType = sum->InputDistType(1);

  if (inType.m_notReped != outType.m_notReped) {
    if (m_dim >= inType.m_numDims)
      return false;
    else
      return sum->m_sumDims.size() > 1;
  }
  else {
    if (m_dim >= outType.m_numDims)
      return false;

    if (sum->m_sumDims.size() <=1)
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
		      return true;
		    }
		  }
		  return false;
		}
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


  if (inType.m_notReped != outType.m_notReped) {
    DistEntry inEntry = inType.m_dists[m_dim];

    //    cout << "applying to sum scatter " << sum->InputDistType(0).PrettyStr() << " -> "
    //	 <<  sum->GetDistType(0).PrettyStr() << endl;


    EntrySet sums = sum->m_sumDims;
    if (!sums.erase(inEntry)) {
      cout << "scalar " << sum->GetNameStr(0) << endl;
      cout << "didn't find entry\n";
      cout << "in " << inType.str()
	   << "\nout " << outType.str()
	   << "\ninEntry: " << inEntry.str() << endl;
      throw;
    }

    if (sums.empty())
      throw;

    DistType intType = sum->GetDistType(0);
    {
      DimSet set = inType.m_notReped.DistEntryDimSet();
      DimSet tempSet = inEntry.DistEntryDimSet();
      set.insert(tempSet.begin(), tempSet.end()); // order them correctly
      DimVec vec;
      vec.insert(vec.end(), set.begin(), set.end());
      intType.m_notReped.DimsToDistEntry(vec);
    }

    //    cout << "intType " << intType.PrettyStr() << endl;

    TempVarNode *temp = new TempVarNode(intType, sums, sum->GetInputName(0).m_name);
    temp->AddInput(node->Input(1), node->InputConnNum(1));
    poss->AddNode(temp);

    //    cout << "tempType " << temp->GetDistType(0).PrettyStr() << endl;

    EntrySet set;
    set.insert(inEntry);
    
    SumScatterUpdateNode *newSum = new SumScatterUpdateNode(COEFONE, set);
    newSum->AddInput(node->Input(0), node->InputConnNum(0));
    newSum->AddInput(temp, 0);
    poss->AddNode(newSum);


    node->ChangeInput2Way(node->Input(0), node->InputConnNum(0),
			  newSum, 0);

    //    cout << "new sum scatter " << newSum->InputDistType(0).PrettyStr() << " -> "
    //	 <<  newSum->GetDistType(0).PrettyStr() << endl;



    if (sums.empty()) {
      cout << "sums is empty1\n";
      throw;
    }
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
	  
	  if (outVec != tempVec)
	    needRedist = true;
	  /*
	  DimVecIter outIter = outVec.begin();
	  DimVecIter tempIter = tempVec.begin();
	  for(; !needRedist && outIter != outVec.end(); ++outIter, ++tempIter) {
	    if (*outIter != *tempIter)
	      needRedist = true;
	  }
	  */

	

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

	  EntrySet sums = sum->m_sumDims;
	  EntrySetIter sumIter = sums.begin();
	  for(; sumIter != sums.end(); ++sumIter) {
	    if (*sumIter == *iter) {
	      sums.erase(sumIter);
	      break;
	    }
	  }


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

	    if (!newOutTypeFull.IsSane()) {
	      cout << "created new " << newSum->GetDistType(0).str() << " -> " << newOutTypeFull.str() << endl;
	      throw;
	    }
	    if (newSum->GetDistType(0).m_numDims != newOutTypeFull.m_numDims)
	      throw;
	  }

	  if (DistTypeNotEqual(newOutTypeFull, outType)) {
	    if (needRedist)
	      node->ChangeInput2Way(node->Input(0), node->InputConnNum(0),
				    newRedist, 0);
	    else
	      node->ChangeInput2Way(node->Input(0), node->InputConnNum(0),
				    newSum, 0);
	    if (sums.empty()) {
	      cout << "sums is empty2\n";
	      cout << newOutTypeFull.str() << " != " << outType.str() << endl;
	      throw;
	    }

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
    throw;
  }
}

bool MoveSumScatterRedistAfter::CanApply(const Poss *poss, const Node *node) const
{
  const SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;
  
  const DistType &inType = sum->InputDistType(0);
  const DistType &outType = sum->InputDistType(1);

  for (Dim dim = 0; dim < outType.m_numDims; ++dim) {
    DistEntry inEntry = inType.m_dists[dim];
    DistEntry outEntry = outType.m_dists[dim];
    
    if (inEntry != outEntry) {
      DimVec inVec = inEntry.DistEntryDims();
      DimSet inSet;
      inSet.insert(inVec.begin(), inVec.end());
      DimVec outVec = outEntry.DistEntryDims();
      DimSet outSet;
      outSet.insert(outVec.begin(), outVec.end());

	if (includes(outSet.begin(), outSet.end(),
		     inSet.begin(), inSet.end())) 
	  {
	    bool diff = false;
	    DimVecIter inIter = inVec.begin();
	    DimVecIter outIter = outVec.begin();
	    for(; !diff && inIter != inVec.end(); ++inIter, ++outIter) {
	      if (*inIter != *outIter)
		diff = true;
	    }
	    
	    if (diff) {
	      EntrySetIter iter = sum->m_sumDims.begin();
	      for(; iter != sum->m_sumDims.end(); ++iter) {
		DistEntry sumEntry = *iter;
		DimVec sumVec = sumEntry.DistEntryDims();
		DimSet sumSet;
		sumSet.insert(sumVec.begin(), sumVec.end());
		if (outSet.size() == (inSet.size() + sumSet.size())
		    && includes(outSet.begin(), outSet.end(),
				sumSet.begin(), sumSet.end())) 
		  {
		    return true;
		  }
	      }
	    }
	  }
    }
  }
  return false;
}

void MoveSumScatterRedistAfter::Apply(Poss *poss, Node *node) const
{
  SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;
  
  const DistType &inType = sum->InputDistType(0);
  const DistType &outType = sum->InputDistType(1);
  DistType outTypeInt = outType;

  for (Dim dim = 0; dim < outTypeInt.m_numDims; ++dim) {
    DistEntry inEntry = inType.m_dists[dim];
    DistEntry outEntry = outType.m_dists[dim];
    
    if (inEntry != outEntry) {
      DimVec inVec = inEntry.DistEntryDims();
      DimSet inSet;
      inSet.insert(inVec.begin(), inVec.end());
      DimVec outVec = outEntry.DistEntryDims();
      DimSet outSet;
      outSet.insert(outVec.begin(), outVec.end());

	if (includes(outSet.begin(), outSet.end(),
		     inSet.begin(), inSet.end())) 
	  {
	    bool diff = false;
	    DimVecIter inIter = inVec.begin();
	    DimVecIter outIter = outVec.begin();
	    for(; !diff && inIter != inVec.end(); ++inIter, ++outIter) {
	      if (*inIter != *outIter)
		diff = true;
	    }
	    
	    if (diff) {
	      EntrySetIter iter = sum->m_sumDims.begin();
	      for(; iter != sum->m_sumDims.end(); ++iter) {
		DistEntry sumEntry = *iter;
		DimVec sumVec = sumEntry.DistEntryDims();
		DimSet sumSet;
		sumSet.insert(sumVec.begin(), sumVec.end());
		if (outSet.size() == (inSet.size() + sumSet.size())
		    && includes(outSet.begin(), outSet.end(),
				sumSet.begin(), sumSet.end())) 
		  {
		    inVec.insert(inVec.end(), sumVec.begin(), sumVec.end());
		    outTypeInt.m_dists[dim].DimsToDistEntry(inVec);
		  }
	      }
	    }
	  }
    }
  }

  TempVarNode *temp = new TempVarNode(outTypeInt, node->GetInputName(0).m_name);
  temp->AddInput(node->Input(1), node->InputConnNum(1));
  poss->AddNode(temp);

  SumScatterUpdateNode *newSum = new SumScatterUpdateNode(sum->m_coef, sum->m_sumDims);
  newSum->AddInput(node->Input(0), node->InputConnNum(0));
  newSum->AddInput(temp, 0);
  poss->AddNode(newSum);


  RedistNode *redist = new RedistNode(outType);
  redist->AddInput(newSum, 0);
  poss->AddNode(redist);


  sum->RedirectChildren(redist, 0);
	    node->m_poss->DeleteChildAndCleanUp(node);
}

#endif




