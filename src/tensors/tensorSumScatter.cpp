/*
  This file is part of DxTer.
  DxTer is a prototype using the Design by Transformation (DxT)
  approach to program generation.
  
  Copyright (C) 2014, The University of Texas and Bryan Marker
  
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
#include "tensorRedist.h"
#include <ostream>
#include <sstream>
#include <algorithm>
#include "helperNodes.h"
#include "tunnel.h"

#if DOTENSORS

void GetMultiModeScatterInfo(const DistType &srcType,
			     const DistType &destType,
			     const EntryList &m_sumDims,
			     DimVec *redDims,
			     DimVec *scatDims)
{
  EntrySet sumDimSet;
  sumDimSet.insert(m_sumDims.begin(),m_sumDims.end());
  
  for (int dim = srcType.m_numDims-1; dim >= 0; --dim) {
    if (sumDimSet.find(srcType.m_dists[dim]) != sumDimSet.end()) {
      redDims->push_back(dim);
    }
  }
  if (redDims->size() != m_sumDims.size())
    throw;
  

  DimVecIter redDimsIter = redDims->begin();
  for(; redDimsIter != redDims->end(); ++redDimsIter) {
    DistEntry entry = srcType.m_dists[*redDimsIter];
    DimVec entryDims = entry.DistEntryDims();
    for(Dim dim = 0; dim < destType.m_numDims; ++dim) {
      DimVec destDims = destType.m_dists[dim].DistEntryDims();
      DimVec srcDims = srcType.m_dists[dim].DistEntryDims();
      srcDims.insert(srcDims.end(), entryDims.begin(), entryDims.end());
      if (destDims == srcDims) {
	scatDims->push_back(dim);
	break;
      }
    }
  }
  if (scatDims->size() != redDims->size())
    throw;



}

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

SumScatterUpdateNode::SumScatterUpdateNode(Coef coef, const EntryList &sumDims)
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
  string tmp = "SumScatterUpdate";
  EntryListConstIter iter = m_sumDims.begin();
  for(; iter != m_sumDims.end(); ++iter) {
    tmp += "|";
    tmp += char(((*iter).m_val)+48);
  }
  return tmp;
}

void SumScatterUpdateNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    m_cost = 0;

    //need to update everything
    if (InputDataType(1).HasPerm())
      throw;

    const DistType &inType = InputDataType(0).GetDist();
    const DistType &outType = InputDataType(1).GetDist();


    
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
      EntryListConstIter iter = m_sumDims.begin();
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
      EntryListIter iter = m_sumDims.begin();
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
    
    unsigned int numProcs = 1;

    EntryListIter iter = m_sumDims.begin();
    for(; iter != m_sumDims.end(); ++iter) {
      DimVec sums = (*(iter)).DistEntryDims();
      DimVecConstIter iter = sums.begin();
      for(; iter != sums.end(); ++iter) {
	numProcs *= GridLens[*iter];
      }
    }

    DLANode *input = (DLANode*)(Input(1));
    ConnNum num = InputConnNum(1);
    const Dim numDims = input->NumDims(num);

    if (inType.m_notReped != outType.m_notReped) {
      input = (DLANode*)(Input(0));
      num = InputConnNum(0);
      m_cost = input->LocalLen(num,0)->NumSizes() * ReduceScatter(numProcs,numProcs);
      m_cost += (PSIR+PSIW)*(numProcs+1);
    }
    else {
      const unsigned int totNumIters = input->LocalLen(num,0)->NumSizes();
      for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
	Cost temp = 1;
	for (Dim dim = 0; dim < numDims; ++dim) {
	  temp *= (*(input->LocalLen(num,dim)))[iteration];
	}
	m_cost += ReduceScatter(temp*numProcs, numProcs);
	m_cost += (PSIR+PSIW)*(temp*numProcs + temp);
      }
    }
  }
}

void SumScatterUpdateNode::CheckSumDimsInOutput() const
{
  const DistType &outType = InputDataType(1).GetDist();
    
  if (InputDataType(0).GetDist().m_notReped != outType.m_notReped)
      return;

    DimSet allSumDims;
    EntryListConstIter iter = m_sumDims.begin();
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
      cout << "inType " << InputDataType(0).GetDist().str() << endl;
      m_poss->PrintTransVec();
      throw;
    }
}

Phase SumScatterUpdateNode::MaxPhase() const
{

  const DistType &inType = InputDataType(0).GetDist();
  const DistType &outType = InputDataType(1).GetDist();

  if (CurrPhase >= SUMSCATTERTENSORPHASE) {
    DistEntry sumDims;
    if (m_sumDims.size() != 1) {
#if !ALLOWMULTIMODESCATTER
      return SUMSCATTERTENSORPHASE;
#else
      if (inType.m_notReped == outType.m_notReped) {
	EntrySet sumSet;
	sumSet.insert(m_sumDims.begin(), m_sumDims.end());

	for (Dim dim = 0; !sumSet.empty() && dim < outType.m_numDims; ++dim) {
	  DistEntry inEntry = inType.m_dists[dim];
	  DistEntry outEntry = outType.m_dists[dim];
	  if (inEntry != outEntry) {
	    if (inEntry.IsStar()) {
	      //[*] -> ...
	      if (!sumSet.erase(outEntry)) {
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
	      DimVec outVec = outEntry.DistEntryDims();
	      
	      if (inVec[0] != outVec[0])
		return SUMSCATTERTENSORPHASE;

	      DimVec suff;
	      GetSuffix(inVec, outVec,
			suff);

	      DistEntry entry;
	      entry.DimsToDistEntry(suff);

	      if (!sumSet.erase(outEntry))
		return SUMSCATTERTENSORPHASE;	
	    }
	  }
	}
	if (!sumSet.empty())
	  throw;
      }
#endif
    }
    else {
      if (inType.m_notReped == outType.m_notReped) {
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
  }

  return NUMPHASES;
}

void SumScatterUpdateNode::PrintCode(IndStream &out)
{
  const DistType &m_srcType = InputDataType(0).GetDist();
  const DistType &m_destType = InputDataType(1).GetDist();

  Dim srcNumDims = m_srcType.m_numDims;

  string inName = GetInputName(0).str();
  string outName = GetName(0).str();

  if (CurrPhase <= SUMSCATTERTENSORPHASE)
    return;

#if !ALLOWMULTIMODESCATTER
  if (m_sumDims.size() != 1)
    throw;
#else
  if (m_sumDims.empty())
    throw;
#endif


  if (m_sumDims.size() == 1) {
    DistEntry sumDims = *(m_sumDims.begin());
    
    out.Indent();
    *out << "   // " << GetName(0).PrettyStr() 
	 << " <- " << GetInputName(0).PrettyStr() 
	 << " (with SumScatter on " << sumDims.PrettyStr() << ")\n";
    
    out.Indent();

    if (m_srcType.m_notReped != m_destType.m_notReped) {
      if (m_coef != COEFZERO) {
	*out << outName << ".ReduceToOneUpdateRedistFrom( "
	     << inName << ", ";
	out << m_coef;
	*out << ", ";
      }
      else {
	*out << outName << ".ReduceToOneRedistFrom( "
	     << inName << ", ";
      }
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

      *out << outName;
      if (m_coef == COEFZERO) 
	*out << ".ReduceScatterRedistFrom( ";
      else
	*out << ".ReduceScatterUpdateRedistFrom( ";

      *out << inName << ", ";
    
      if (m_coef != COEFZERO) {
	out << m_coef;
	*out << ", ";
      }

      DimVec scatDims = sumDims.DistEntryDims();


      Dim redDim = srcNumDims;
      for (int dim = srcNumDims-1; dim >= 0; --dim) {
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
  else {
    out.Indent();
    *out << "   // " << GetName(0).PrettyStr() 
	 << " <- " << GetInputName(0).PrettyStr() 
	 << " (with SumScatter on ";
    EntryListIter iter = m_sumDims.begin();
    for(; iter != m_sumDims.end(); ++iter) {
      *out << "(" << (*iter).PrettyStr() << ")";
    }
    *out << ")\n";
    
    out.Indent();

    if (m_srcType.m_notReped != m_destType.m_notReped) {
      if (m_coef != COEFZERO) {
	*out << outName << ".ReduceToOneUpdateRedistFrom( "
	     << inName << ", ";
	out << m_coef;
	*out << ", ";
      }
      else {
	*out << outName << ".ReduceToOneRedistFrom( "
	     << inName << ", ";
      }

      EntrySet dimSet;
      dimSet.insert(m_sumDims.begin(), m_sumDims.end());
      DimVec dims;

      for(Dim dim = 0; dim < srcNumDims; ++dim) {
	if (dimSet.erase(m_srcType.m_dists[dim])) {
	  dims.push_back(dim);
	}
      }
      if (!dimSet.empty())
	throw;


      *out << ModeArrayVarName(dims) << " );\n";
    }
    else {
      if (srcNumDims != (m_destType.m_numDims+m_sumDims.size()))
	throw;

      *out << outName;
      if (m_coef == COEFZERO) 
	*out << ".ReduceScatterRedistFrom( ";
      else
	*out << ".ReduceScatterUpdateRedistFrom( ";

      *out << inName << ", ";
    
      if (m_coef != COEFZERO) {
	out << m_coef;
	*out << ", ";
      }

      DimVec redDims, scatDims;

      GetMultiModeScatterInfo(m_srcType, m_destType,
			      m_sumDims,
			      &redDims,
			      &scatDims);
	

      *out << ModeArrayVarName(redDims) << ", " << ModeArrayVarName(scatDims) << " );\n";
    }
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

void SumScatterUpdateNode::AddVariables(VarSet &set) const
{
  DLAOp<2,1>::AddVariables(set);
#if ALLOWMULTIMODESCATTER
  if (m_sumDims.size() > 1) {
    const DistType &m_srcType = InputDataType(0).GetDist();
    const DistType &m_destType = InputDataType(1).GetDist();
    if (m_srcType.m_notReped != m_destType.m_notReped) {
      EntrySet dimSet;
      dimSet.insert(m_sumDims.begin(), m_sumDims.end());
      DimVec dims;

      for(Dim dim = 0; dim < m_srcType.m_numDims; ++dim) {
	if (dimSet.erase(m_srcType.m_dists[dim])) {
	  dims.push_back(dim);
	}
      }
      if (!dimSet.empty())
	throw;
      Var var(ModeArrayVarType, dims);
      set.insert(var);
    }
    else {
      DimVec redDims, scatDims;

      GetMultiModeScatterInfo(m_srcType, m_destType,
			      m_sumDims,
			      &redDims,
			      &scatDims);
      {
	Var var(ModeArrayVarType, redDims);
	set.insert(var);
      }
      {
	Var var(ModeArrayVarType, scatDims);
	set.insert(var);
      }
    }
  }
#else
  if (m_sumDims.size() > 1) 
    throw;
#endif 
}



bool SeparateRedistFromSumScatter::CanApply(const Node *node) const
{
  const SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;
  const DistType &inType = sum->InputDataType(0).GetDist();
  const DistType &outType = sum->InputDataType(1).GetDist();

  if (inType.m_notReped != outType.m_notReped)
    return false;
  

  const EntryList &sumDims = sum->m_sumDims;

  DimSet allSumDims;
  EntryListConstIter iter = sumDims.begin();
  for(; iter != sumDims.end(); ++iter) {
    DimVec vec = (*iter).DistEntryDims();
    allSumDims.insert(vec.begin(), vec.end());
  }

  for(Dim dim = 0; dim < inType.m_numDims; ++dim) {
    DistEntry inEntry = inType.m_dists[dim];
    if (dim >= outType.m_numDims) {
      cout << "dim >= outType.m_numDims\n";
      cout << sum->InputDataType(0).GetDist().str() << " -> " 
	   << outType.str() << " with " << sum->m_sumDims.size() << endl;
      sum->m_poss->PrintTransVec();
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

void SeparateRedistFromSumScatter::Apply(Node *node) const
{
  SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;

  DistType inTypeInt = sum->InputDataType(0).GetDist();
  const DistType &outType = sum->InputDataType(1).GetDist();

  const EntryList &sumDims = sum->m_sumDims;

  DimSet allSumDims;
  EntryListConstIter iter = sumDims.begin();
  for(; iter != sumDims.end(); ++iter) {
    DimVec vec = (*iter).DistEntryDims();
    allSumDims.insert(vec.begin(), vec.end());
  }

  for(Dim dim = 0; dim < inTypeInt.m_numDims; ++dim) {
    DistEntry inEntry = inTypeInt.m_dists[dim];
    if (dim >= outType.m_numDims) {
      cout << "dim >= outType.m_numDims\n";
      cout << sum->InputDataType(0).GetDist().str() << " -> " 
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

  if (sum->Input(1)->m_poss != node->m_poss) {
    cout << "input 1 on " << sum->Input(1)->m_poss << " and this is on " << node->m_poss << endl;
    throw;
  }
  if (sum->m_poss != node->m_poss) {
    cout << "one\n";
    throw;
  }
  
  if (DistTypeEqual(inTypeInt,sum->InputDataType(0).GetDist())) {
    cout << "same dist\n";
    throw;
  }


  RedistNode *redist = new RedistNode(inTypeInt);
  redist->AddInput(sum->Input(0), sum->InputConnNum(0));

  if (inTypeInt == sum->InputDataType(0).GetDist())
    throw;
  
  node->m_poss->AddNode(redist);

  if (!inTypeInt.IsSane()) {
    cout << "!sane: " << inTypeInt.str() << endl;
    cout << "was " << sum->InputDataType(0).GetDist().str() << " -> "
	 << outType.str() << " with " << sumDims.size() << endl;
      EntryListConstIter iter = sumDims.begin();
      for(; iter != sumDims.end(); ++iter) {
	DimVec vec = (*iter).DistEntryDims();
	cout << "sum on " << *(vec.begin()) << endl;
      }

    throw;
  }
  
  SumScatterUpdateNode *sum2 = new SumScatterUpdateNode(sum->m_coef, sumDims);
  sum2->AddInput(redist, 0);
  sum2->AddInput(sum->Input(1), sum->InputConnNum(1));
  node->m_poss->AddNode(sum2);


  sum->RedirectChildren(sum2,0);

  node->m_poss->DeleteChildAndCleanUp(node);
}

bool SplitSumScatter::CanApply(const Node *node) const
{
  const SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;
  
  const DistType &inType = sum->InputDataType(0).GetDist();
  const DistType &outType = sum->InputDataType(1).GetDist();

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

    const EntryList &sumDims = sum->m_sumDims;
  


    DimVec inVec = inEntry.DistEntryDims();
    DimSet inSet;
    inSet.insert(inVec.begin(), inVec.end());
    DimVec outVec = outEntry.DistEntryDims();
    DimSet outSet;
    outSet.insert(outVec.begin(), outVec.end());

	
    bool found = false;
    EntryListConstIter iter = sumDims.begin();
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

void SplitSumScatter::Apply(Node *node) const
{
  SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;

  const DistType &inType = sum->InputDataType(0).GetDist();
  const DistType &outType = sum->InputDataType(1).GetDist();


  if (inType.m_notReped != outType.m_notReped) {
    DistEntry inEntry = inType.m_dists[m_dim];

    //    cout << "applying to sum scatter " << sum->InputDataType(0).m_dist.PrettyStr() << " -> "
    //	 <<  sum->GetDistType(0).PrettyStr() << endl;

    EntryList sums = sum->m_sumDims;
    EntryListIter iter = sums.begin();
    bool found = false;
    for(; iter != sums.end() && !found; ++iter) {
      if (*iter == inEntry) {
	sums.erase(iter);
	found = true;
	break;
      }
    }
    if (!found)
      throw;


    DistType intType = sum->DataType(0).GetDist();
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
    node->m_poss->AddNode(temp);

    //    cout << "tempType " << temp->GetDistType(0).PrettyStr() << endl;

    EntryList list;
    list.push_back(inEntry);
    
    SumScatterUpdateNode *newSum = new SumScatterUpdateNode(COEFZERO, list);
    newSum->AddInput(node->Input(0), node->InputConnNum(0));
    newSum->AddInput(temp, 0);
    node->m_poss->AddNode(newSum);


    node->ChangeInput2Way(node->Input(0), node->InputConnNum(0),
			  newSum, 0);

    //    cout << "new sum scatter " << newSum->InputDataType(0).GetDist().PrettyStr() << " -> "
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
	
    EntryListIter iter = sum->m_sumDims.begin();
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

	  EntryList sums = sum->m_sumDims;
	  EntryListIter sumIter = sums.begin();
	  for(; sumIter != sums.end(); ++sumIter) {
	    if (*sumIter == *iter) {
	      sums.erase(sumIter);
	      break;
	    }
	  }


	  TempVarNode *temp = new TempVarNode(newOutTypeNoSums, sums, node->GetInputName(0).m_name);
	  temp->AddInput(node->Input(1), node->InputConnNum(1));
	  node->m_poss->AddNode(temp);
	  

	  EntryList list;
	  list.push_back(*iter);

	  SumScatterUpdateNode *newSum = new SumScatterUpdateNode(COEFZERO, list);
	  newSum->AddInput(node->Input(0), node->InputConnNum(0));
	  newSum->AddInput(temp, 0);
	  node->m_poss->AddNode(newSum);

	  RedistNode *newRedist = NULL;
	  if (needRedist) {
	    if (newOutTypeFull == newSum->DataType(0).GetDist())
	      needRedist = false;
	    else {
	      newRedist = new RedistNode(newOutTypeFull);
	      newRedist->AddInput(newSum, 0);
	      
	      node->m_poss->AddNode(newRedist);
	      
	      if (!newOutTypeFull.IsSane()) {
		cout << "created new " << newSum->DataType(0).GetDist().str() << " -> " << newOutTypeFull.str() << endl;
		throw;
	      }
	      if (newSum->DataType(0).GetDist().m_numDims != newOutTypeFull.m_numDims)
		throw;
	    }
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

bool MoveSumScatterRedistAfter::CanApply(const Node *node) const
{
  const SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;
  
  const DistType &inType = sum->InputDataType(0).GetDist();
  const DistType &outType = sum->InputDataType(1).GetDist();

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
	      EntryListConstIter iter = sum->m_sumDims.begin();
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

void MoveSumScatterRedistAfter::Apply(Node *node) const
{
  SumScatterUpdateNode *sum = (SumScatterUpdateNode*)node;
  
  const DistType &inType = sum->InputDataType(0).GetDist();
  const DistType &outType = sum->InputDataType(1).GetDist();
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
	      EntryListIter iter = sum->m_sumDims.begin();
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

  RedistNode *redist = new RedistNode(outTypeInt);
  redist->AddInput(node->Input(1), node->InputConnNum(1));
  node->m_poss->AddNode(redist);

  if (outTypeInt == redist->InputDataType(0).GetDist())
    throw;

  SumScatterUpdateNode *newSum = new SumScatterUpdateNode(sum->m_coef, sum->m_sumDims);
  newSum->AddInput(node->Input(0), node->InputConnNum(0));
  newSum->AddInput(redist, 0);
  node->m_poss->AddNode(newSum);


  RedistNode *redist2 = new RedistNode(outType);
  redist2->AddInput(newSum, 0);
  node->m_poss->AddNode(redist2);

  if (outType == redist2->InputDataType(0).GetDist())
    throw;
  
  
  sum->RedirectChildren(redist2, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#endif

