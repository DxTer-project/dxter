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

#if DOTENSORS

#include "tensorSumScatter.h"


DimVec GetCommonPrefix(const DimVec &dims1, const DimVec &dims2);
void GetCommonPrefix(const DimVec &dims1, const DimVec &dims2, 
		     DimVec &pref,
		     DimVec &suff1, DimVec &suff2);
void GetCommonSuffix(const DimVec &dims1, const DimVec &dims2, 
		     DimVec &suff,
		     DimVec &pref1, DimVec &pref2);
void GetDifferentSuffix(const DimVec &dims1, const DimVec &dims2, 
			DimVec &suff1, DimVec &suff2);
void GetSuffix(const DimVec &dims1, const DimVec &dims2, 
	       DimVec &suff);

RedistNode::RedistNode() 
{
  m_info.m_dist.SetToDefault(0);
  m_lsizes = NULL;
}

RedistNode::RedistNode(const DistType &destType)
{
  m_info.m_dist = destType;
  m_lsizes = NULL;
}

RedistNode::~RedistNode()
{
  if (m_lsizes) {
    if (m_isArray)
      delete [] m_lsizes;
    else
      delete m_lsizes;
    m_lsizes = NULL;
  }
}

void RedistNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig,shallow, possMerging);
  const RedistNode *origNode = (RedistNode*)orig;
  m_info = origNode->m_info;
}

NodeType RedistNode::GetType() const
{
  if (m_name.length()) 
    return m_name;
  else {
    return (string)"RedistNode to " +  m_info.m_dist.QuickStr();
  }
}

void RedistNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLANode::Prop();

    if (m_inputs.size() != 1) {
      cout << "m_inputs.size() != 1\n";
      throw;
    }
    
    if (!m_children.size())
      throw;
  
    if (!m_name.length())
      m_name = (string)"RedistNode to " +  m_info.m_dist.QuickStr();
    DLANode *parent = (DLANode*)Input(0);
    const DistType &m_srcType = parent->DataType(InputConnNum(0)).m_dist;
    parent->Prop();
    const Dim numDims = m_info.m_dist.m_numDims;
    
    if (numDims != InputNumDims(0)) {
      cout << "numDims " << numDims << " " << InputNumDims(0) << endl;
      cout << "input " << Input(0)->GetNodeClass() << endl;
      cout << GetInputNameStr(0) << endl;
      cout << m_info.m_dist.PrettyStr() << endl;
      throw;
    }
    if (m_srcType.m_numDims != numDims)
      throw;

    if (m_srcType.m_numDims != numDims) {
      cout << m_info.m_dist.str() << " <- " << m_srcType.str() << endl;
      throw;
    }

    if (!m_info.m_dist.IsSane()) {
      cout << m_info.m_dist.str() << endl;
      m_poss->PrintTransVec();
      throw;
    }

    DimSet diffs;
    
    for (Dim dim = 0; dim < numDims; ++dim) {
      if (m_srcType.m_dists[dim] != m_info.m_dist.m_dists[dim]) {
	diffs.insert(dim);
      }
    }

    if (diffs.empty()) {
      //      throw;
      m_cost = 0;
    }
    else if (diffs.size() == 1) {
      const Dim dim = *(diffs.begin());
      DimVec src = m_srcType.m_dists[dim].DistEntryDims();
      DimVec dest = m_info.m_dist.m_dists[dim].DistEntryDims();

      
      if (src.empty() || IsPrefix(src, dest)) {
	//local memory copy
	//	throw;
	m_cost = 0;
      }
      else if (IsPrefix(dest, src) || dest.empty()) {
	m_cost = 0;
	const unsigned int totNumIters = m_lsizes[0].NumSizes();
	unsigned int numProcs = 1;
	DimVecIter iter = src.begin() + dest.size();
	for(; iter != src.end(); ++iter) {
	  numProcs *= GridLens[*iter];
	}
	for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
	  Cost temp = 1;
	  for (Dim dim = 0; dim < numDims; ++dim) {
	    temp *= m_lsizes[dim][iteration];
	  }
	  m_cost += AllGather(temp * numProcs, numProcs);
	}
      }
      else {
	m_cost = 0;
	const unsigned int totNumIters = m_lsizes[0].NumSizes();
	unsigned int numProcs = 1;
	DimVecIter iter = src.begin();
	DimSet unionSet;
	for(; iter != src.end(); ++iter) {
	  numProcs *= GridLens[*iter];
	  unionSet.insert(*iter);
	}
	iter = dest.begin();
	for(; iter != dest.end(); ++iter) {
	  if (unionSet.find(*iter) == unionSet.end())
	    numProcs *= GridLens[*iter];
	}

	for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
	  Cost temp = 1;
	  for (Dim dim = 0; dim < numDims; ++dim) {
	    temp *= m_lsizes[dim][iteration];
	  }
	  m_cost += AllToAll(temp * numProcs, numProcs);
	}
      }
    }
    else if (diffs.size() > 2) {
      m_cost = 0;
    }
    else {

      DimSetConstIter iterDiffs = diffs.begin();
      Dim diff1 = *iterDiffs;
      ++iterDiffs;
      Dim diff2 = *iterDiffs;

      DimVec srcDims = m_srcType.m_dists[diff1].DistEntryDims();
      DimVec destDims = m_info.m_dist.m_dists[diff1].DistEntryDims();
      DimVec otherSrcDims = m_srcType.m_dists[diff2].DistEntryDims();
      DimVec otherDestDims = m_info.m_dist.m_dists[diff2].DistEntryDims();
    
      DimVec commonSuff1, commonSuff2, otherSrcPref, destPref, srcPref, otherDestPref;
      GetCommonSuffix(otherSrcDims, destDims, commonSuff1, otherSrcPref, destPref);
      GetCommonSuffix(srcDims, otherDestDims, commonSuff2, srcPref, otherDestPref);

    
      m_cost = 0;
      const unsigned int totNumIters = m_lsizes[0].NumSizes();
      unsigned int numProcs = 1;
      DimSet unionSet;
      
      DimVecIter iter = commonSuff1.begin();
      for(; iter != commonSuff1.end(); ++iter) {
	if (unionSet.insert(*iter).second)
	  numProcs *= GridLens[*iter];
      }
      iter = commonSuff2.begin();
      for(; iter != commonSuff2.end(); ++iter) {
	if (unionSet.insert(*iter).second)
	  numProcs *= GridLens[*iter];
      }

      for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
	Cost temp = 1;
	for (Dim dim = 0; dim < numDims; ++dim) {
	  temp *= m_lsizes[dim][iteration];
	}
	m_cost += AllToAll(temp * numProcs, numProcs);
      }
    }
  }
}


DimVec GetCommonPrefix(const DimVec &dims1, const DimVec &dims2)
{
  DimVec pref;
  DimVecConstIter iter1 = dims1.begin();
  DimVecConstIter iter2 = dims2.begin();
  while(iter1 != dims1.end() && iter2 != dims2.end()) {
    if (*iter1 == *iter2) {
      pref.push_back(*iter1);
      ++iter1;
      ++iter2;
    }
    else {
      break;
    }
  }
  return pref;
}

void GetCommonSuffix(const DimVec &dims1, const DimVec &dims2, 
		     DimVec &suff,
		     DimVec &pref1, DimVec &pref2)
{
  DimVecConstRevIter iter1 = dims1.rbegin();
  DimVecConstRevIter iter2 = dims2.rbegin();
  while(iter1 != dims1.rend() && iter2 != dims2.rend()) {
    if (*iter1 == *iter2) {
      suff.insert(suff.begin(), *iter1);
      ++iter1;
      ++iter2;
    }
    else {
      break;
    }
  }
  while (iter1 != dims1.rend()) {
    pref1.insert(pref1.begin(),*iter1);
    ++iter1;
  }
  while (iter2 != dims2.rend()) {
    pref2.insert(pref2.begin(),*iter2);
    ++iter2;
  }
}

void GetCommonPrefix(const DimVec &dims1, const DimVec &dims2, 
		     DimVec &pref,
		     DimVec &suff1, DimVec &suff2)
{
  DimVecConstIter iter1 = dims1.begin();
  DimVecConstIter iter2 = dims2.begin();
  while(iter1 != dims1.end() && iter2 != dims2.end()) {
    if (*iter1 == *iter2) {
      pref.push_back(*iter1);
      ++iter1;
      ++iter2;
    }
    else {
      break;
    }
  }
  suff1.insert(suff1.begin(), iter1, dims1.end());
  suff2.insert(suff2.begin(), iter2, dims2.end());
}

void GetDifferentSuffix(const DimVec &dims1, const DimVec &dims2, 
			DimVec &suff1, DimVec &suff2)
{
  DimVecConstIter iter1 = dims1.begin();
  DimVecConstIter iter2 = dims2.begin();
  while (iter1 != dims1.end() && iter2 != dims2.end()) {
    if (*iter1 == *iter2) {
      ++iter1;
      ++iter2;
    }
    else
      break;
  }
  if (iter1 == dims1.end())
    throw;
  if (iter2 == dims2.end())
    throw;
  while (iter1 != dims1.end()) {
    suff1.push_back(*iter1);
    ++iter1;
  }
  while (iter2 != dims2.end()) {
    suff2.push_back(*iter2);
    ++iter2;
  }
}

void GetSuffix(const DimVec &dims1, const DimVec &dims2, 
			DimVec &suff)
{
  DimVecConstIter iter1 = dims1.begin();
  DimVecConstIter iter2 = dims2.begin();
  while (iter1 != dims1.end() && iter2 != dims2.end()) {
    if (*iter1 == *iter2) {
      ++iter1;
      ++iter2;
    }
    else
      break;
  }
  if (iter2 == dims2.end())
    throw;
  while (iter2 != dims2.end()) {
    suff.push_back(*iter2);
    ++iter2;
  }
}


Phase RedistNode::MaxPhase() const
{
  DLANode *parent = (DLANode*)Input(0);
  const DistType &m_srcType = parent->DataType(InputConnNum(0)).m_dist;
  const Dim numDims = m_info.m_dist.m_numDims;
  
  bool foundMemCopy = false;
  bool foundAllGather = false;
  bool foundAllToAll = false;
  DimSet diffs;
  
  for (Dim dim = 0; dim < numDims; ++dim) {
    DistEntry srcDistEntry = m_srcType.m_dists[dim];
    DistEntry destDistEntry = m_info.m_dist.m_dists[dim];
    if (srcDistEntry != destDistEntry) {
      DimVec srcDims = srcDistEntry.DistEntryDims();
      DimVec destDims = destDistEntry.DistEntryDims();
      if (diffs.size() == 1) {
	foundAllToAll = true;
	diffs.insert(dim);
	Dim otherDim = *(diffs.begin());
	DimVec otherSrcDims = m_srcType.m_dists[otherDim].DistEntryDims();
	DimVec otherDestDims = m_info.m_dist.m_dists[otherDim].DistEntryDims();
	    
	DimVec commonSuff1, commonSuff2, otherSrcPref, destPref, srcPref, otherDestPref;
	GetCommonSuffix(otherSrcDims, destDims, commonSuff1, otherSrcPref, destPref);
	GetCommonSuffix(srcDims, otherDestDims, commonSuff2, srcPref, otherDestPref);
	    
	if (commonSuff1.empty() && commonSuff2.empty())
	  return ROTENSORPHASE;

	if (otherSrcPref != otherDestPref) {
	  return ROTENSORPHASE;
	}
	if (srcPref != destPref) {
	  return ROTENSORPHASE;
	}
      }
      else if (srcDistEntry.IsStar() || IsPrefix(srcDims, destDims)) {
	foundMemCopy = true;
	diffs.insert(dim);
	if (foundAllGather || foundAllToAll) {
	  return ROTENSORPHASE;
	}
      }
      else if (destDistEntry.IsStar() || IsPrefix(destDims, srcDims)) {
	diffs.insert(dim);
	if (foundAllGather || foundMemCopy || foundAllToAll)
	  return ROTENSORPHASE;
	else
	  foundAllGather = true;
      }
      else {
	foundAllToAll = true;
	diffs.insert(dim);
      }
    }
  }

  if (diffs.size() == 1 && foundAllToAll) {
    Dim dim = *(diffs.begin());
    DistEntry srcDistEntry = m_srcType.m_dists[dim];
    DistEntry destDistEntry = m_info.m_dist.m_dists[dim];

    DimSet srcSet = srcDistEntry.DistEntryDimSet();
    DimSet destSet = destDistEntry.DistEntryDimSet();
    if (srcSet.size() != destSet.size())
      return ROTENSORPHASE;
    if (!includes(srcSet.begin(), srcSet.end(),
		  destSet.begin(), destSet.end())) {
      return ROTENSORPHASE;
    }
  }

  return NUMPHASES;
}

const Dim RedistNode::NumDims(unsigned int num) const
{
  if (num > 0)
    throw;
  return InputNumDims(0);
}

const Sizes* RedistNode::Len(unsigned int num, Dim dim) const
{
  if (num > 0)
    throw;
  return InputLen(0,dim);
}

const Sizes* RedistNode::LocalLen(unsigned int num, Dim dim) const
{
  if (num > 0)
    throw;
  if (!m_isArray) {
    return m_lsizes;
  }
  else
    return m_lsizes+dim;
}

void RedistNode::ClearDataTypeCache()
{
  if (m_lsizes) {
    if (!m_isArray) {
      delete m_lsizes;
    }
    else {
      delete [] m_lsizes;
    }
    m_lsizes = NULL;
  }
}

void RedistNode::BuildDataTypeCache()
{
  if (m_lsizes)
    return;

  DLANode *in = (DLANode*)Input(0);
  unsigned int num = InputConnNum(0);
  Dim numDims = in->NumDims(num);
  if (numDims) {
    m_isArray = true;
    m_lsizes = new Sizes[numDims];
    for (Dim dim = 0; dim < numDims; ++dim)
      GetLocalSizes(m_info.m_dist, dim, in->Len(num,dim), m_lsizes+dim);
  }
  else {
    m_isArray = false;
    m_lsizes = new Sizes;
    *m_lsizes = *(in->Len(num,0));
  }
}

void RedistNode::FlattenCore(ofstream &out) const
{
  throw;
  //need to fix flatten code for DistType throughout code base; don't just do a pointer flatten
}
   
void RedistNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  throw;
}

Name RedistNode::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  Name name = GetInputName(0);
  name.m_type = m_info.m_dist;
  return name;
}

void RedistNode::PrintCode(IndStream &out)
{  
  //Reflect in AddVars
  out.Indent();

  const DistType &m_srcType = InputDataType(0).m_dist;
  const Dim numDims = m_info.m_dist.m_numDims;

  string inName = GetInputName(0).str();
  string outName = GetName(0).str(); 

  *out << "   // " << GetName(0).PrettyStr() 
       << " <- " << GetInputName(0).PrettyStr() << endl;
  out.Indent();
  
  if (m_srcType.m_numDims != numDims)
    throw;

  DimSet diffs;

  for (Dim dim = 0; dim < numDims; ++dim) {
    if (m_srcType.m_dists[dim] != m_info.m_dist.m_dists[dim]) {
      diffs.insert(dim);
    }
  }

  if (diffs.empty()) {
    throw;
  }
  else if (diffs.size() == 1) {
    const Dim dim = *(diffs.begin());
    DimVec src = m_srcType.m_dists[dim].DistEntryDims();
    DimVec dest = m_info.m_dist.m_dists[dim].DistEntryDims();
       
    if (src.empty() || IsPrefix(src, dest)) {
      DimVec suff;
      GetSuffix(src, dest, suff);
      /*
       *out << "LocalRedist( " << outName << ", "
       << inName << ", " << dim 
       << ", " << ModeArrayVarName(suff)
       << " );\n";
      */
      *out << outName << ".LocalRedistFrom( " 
	   << inName << ", " << dim 
	   << ", " << ModeArrayVarName(suff)
	   << " );\n";
    }
    else if (IsPrefix(dest, src) || dest.empty()) {
      if (!src.empty()) {
	
	//	*out << "AllGatherRedist( " << outName << ", "
	//	     << inName << ", " << dim << ", ";
	//	DimVec suff;
	//	GetSuffix(dest, src, suff);
	//	*out << ModeArrayVarName(suff) << " );\n";
	*out << outName << ".AllGatherRedistFrom( " << inName
	     << ", " << dim << ", ";
	DimVec suff;
	GetSuffix(dest, src, suff);
	*out << ModeArrayVarName(suff) << " );\n";
      }
      else {
	throw;
	//*out << "AllGatherRedist( " << outName << ", "
	//	     << inName << ", " << dim << " );\n";
      }
    }
    else {
      //      *out << "PermutationRedist( " << outName << ", "
      //      	   << inName << ", " << dim << " );\n";
      *out << outName << ".PermutationRedistFrom( "
      	   << inName << ", " << dim << ", "
	   << ModeArrayVarName(src) << " );\n";
    }
  }
  else if (diffs.size() > 2) {
    throw;
  }
  else {
    DimSetConstIter iter = diffs.begin();
    Dim diff1 = *iter;
    ++iter;
    Dim diff2 = *iter;

    /*
    *out << "AllToAllDoubleIndexRedist( "  << outName << ", "
	 << inName << ", " 
	 << IndexPairVarName(diff1, diff2)
	 << ", ";
    */
    *out << outName << ".AllToAllDoubleIndexRedistFrom( "
	 << inName << ", " 
	 << IndexPairVarName(diff1, diff2)
	 << ", ";


    DimVec firstVec, secondVec;

    DimVec srcDims = m_srcType.m_dists[diff1].DistEntryDims();
    DimVec destDims = m_info.m_dist.m_dists[diff1].DistEntryDims();
    DimVec otherSrcDims = m_srcType.m_dists[diff2].DistEntryDims();
    DimVec otherDestDims = m_info.m_dist.m_dists[diff2].DistEntryDims();
    
    DimVec commonSuff1, commonSuff2, otherSrcPref, destPref, srcPref, otherDestPref;
    GetCommonSuffix(otherSrcDims, destDims, commonSuff1, otherSrcPref, destPref);
    GetCommonSuffix(srcDims, otherDestDims, commonSuff2, srcPref, otherDestPref);

    
    *out << ModeArrayPairVarName(commonSuff2, commonSuff1)
	 << " );\n";
  }
}

void RedistNode::AddVariables(VarSet &set) const
{
  //Reflect in PrintCode
  DLANode::AddVariables(set);
  
  const DistType &m_srcType = InputDataType(0).m_dist;
  const Dim numDims = m_info.m_dist.m_numDims;

  string inName = GetInputName(0).str();
  string outName = GetName(0).str(); 

  if (m_srcType.m_numDims != numDims)
    throw;

  DimSet diffs;
    
  for (Dim dim = 0; dim < numDims; ++dim) {
    if (m_srcType.m_dists[dim] != m_info.m_dist.m_dists[dim]) {
      diffs.insert(dim);
    }
  }

  if (diffs.empty()) {
    throw;
  }
  else if (diffs.size() == 1) {
    const Dim dim = *(diffs.begin());
    DimVec src = m_srcType.m_dists[dim].DistEntryDims();
    DimVec dest = m_info.m_dist.m_dists[dim].DistEntryDims();
       
    if (src.empty() || IsPrefix(src, dest)) {
      DimVec suff;
      GetSuffix(src, dest, suff);
      Var var(suff);
      set.insert(var);
    }
    else if (IsPrefix(dest, src) || dest.empty()) {
      if (!src.empty()) {
	DimVec suff;
	GetSuffix(dest, src, suff);
	Var var(suff);
	set.insert(var);
      }
    }
    else {
      Var var(src);
      set.insert(var);
    }

  }
  else if (diffs.size() > 2) {
    throw;
  }
  else {
    DimSetConstIter iter = diffs.begin();
    Dim diff1 = *iter;
    ++iter;
    Dim diff2 = *iter;

    Var pairVar(diff1, diff2);
    set.insert(pairVar);

    DimVec firstVec, secondVec;

    DimVec srcDims = m_srcType.m_dists[diff1].DistEntryDims();
    DimVec destDims = m_info.m_dist.m_dists[diff1].DistEntryDims();
    DimVec otherSrcDims = m_srcType.m_dists[diff2].DistEntryDims();
    DimVec otherDestDims = m_info.m_dist.m_dists[diff2].DistEntryDims();
    
    DimVec commonSuff1, commonSuff2, otherSrcPref, destPref, srcPref, otherDestPref;
    GetCommonSuffix(otherSrcDims, destDims, commonSuff1, otherSrcPref, destPref);
    GetCommonSuffix(srcDims, otherDestDims, commonSuff2, srcPref, otherDestPref);

    Var arrayPairVar(commonSuff2, commonSuff1);
    set.insert(arrayPairVar);
  }
}


AllReduceNode::AllReduceNode(const DimVec &sumDims, const string &sumIndices)
  : DLAOp<1,1>()
{
  m_sumDims = sumDims;
  m_sumIndices = sumIndices;
}


void AllReduceNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const AllReduceNode *node = (AllReduceNode*)orig;
  DLAOp<1,1>::Duplicate(node, shallow, possMerging);
  m_sumDims = node->m_sumDims;  
  m_sumIndices = node->m_sumIndices;
}

NodeType AllReduceNode::GetType() const
{
  stringstream str;
  str << "AllReduce";
  DimVecConstIter iter = m_sumDims.begin();
  for(; iter != m_sumDims.end(); ++iter)
    str << *iter << ",";
  //  str << m_sumIndices;
  return str.str();
}


void AllReduceNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<1,1>::Prop();


    m_cost = 0;
    unsigned int numProcs = 1;

    DimVecIter iter = m_sumDims.begin();
    for(; iter != m_sumDims.end(); ++iter) {
      numProcs *= GridLens[*iter];
    }
      
    DLANode *input = (DLANode*)(Input(0));
    unsigned int num = InputConnNum(0);

    const unsigned int totNumIters = input->LocalLen(num,0)->NumSizes();
    const Dim numDims = input->NumDims(num);

    for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
      Cost temp = 1;
      for (Dim dim = 0; dim < numDims; ++dim) {
	temp *= (*(input->LocalLen(num,dim)))[iteration];
      }
      m_cost += AllReduce(temp, numProcs);
    }
  }
}

void AllReduceNode::PrintCode(IndStream &out)
{
  throw;
  out.Indent();
  *out << "AllReduceOnDims( " << GetName(0).str()
       << ", m";
  DimVecConstIter iter = m_sumDims.begin();
  for(; iter != m_sumDims.end(); ++iter) {
    *out << "_" << *iter;
  }
  *out << ", " << m_sumIndices << " );\n";
}

void AllReduceNode::FlattenCore(ofstream &out) const
{
  throw;
}

void AllReduceNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  throw;
}



bool RemoveWastedRedist::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  RedistNode *redistNode = (RedistNode*)node;
  const DistType* type = &(redistNode->m_info.m_dist);
  if (node->m_children.size() == 0)
    throw;
  while(redistNode->Input(0) 
        && (redistNode->Input(0)->GetNodeClass() == RedistNode::GetClass()))
    {
      redistNode = (RedistNode*)redistNode->Input(0);
      if (DistTypeEqual(redistNode->m_info.m_dist,*type))
	return true;
      for(unsigned int i = 0; i < redistNode->m_children.size(); ++i) {
	Node *tmp = redistNode->Child(i);
	if (tmp != node && tmp->GetNodeClass() == RedistNode::GetClass()) {
	  if (DistTypeEqual(((RedistNode*)tmp)->m_info.m_dist, *type))
	    return true;
	}
      }
    }
  if (redistNode->Input(0)
      && DistTypeEqual(((DLANode*)(redistNode->Input(0)))->DataType(redistNode->InputConnNum(0)).m_dist,*type))
    return true;
  return false;
}

void RemoveWastedRedist::Apply(Node *node) const
{
  RedistNode *redistNode = (RedistNode*)node;
  const DistType* type = &(redistNode->m_info.m_dist);
  while(redistNode->Input(0) 
        && (redistNode->Input(0)->GetNodeClass() == RedistNode::GetClass()))
    {
      redistNode = (RedistNode*)redistNode->Input(0);
      if (DistTypeEqual(redistNode->m_info.m_dist, *type)) {
	node->RedirectChildren(redistNode, 0);
	node->m_poss->DeleteChildAndCleanUp(node);
	return;
      }
      for(unsigned int i = 0; i < redistNode->m_children.size(); ++i) {
	Node *tmp = redistNode->Child(i);
	if (tmp != node && tmp->GetNodeClass() == RedistNode::GetClass()) {
	  if (DistTypeEqual(((RedistNode*)tmp)->m_info.m_dist,*type)) {
	    node->RedirectChildren(tmp, 0);
	    node->m_poss->DeleteChildAndCleanUp(node);
	    return;
	  }
	}
      }
    }
  if (redistNode->Input(0)
      && DistTypeEqual(((DLANode*)redistNode->Input(0))->DataType(redistNode->InputConnNum(0)).m_dist,*type))
    {
      node->RedirectChildren(redistNode->Input(0), redistNode->InputConnNum(0));
      node->m_poss->DeleteChildAndCleanUp(node);
      return;
    }
  throw;
}


bool RemoveNOPRedistribs::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const DLANode *ddla = (DLANode*)node;
  const Node *parent = ddla->Input(0);
  if (!parent
      || (DistTypeNotEqual(((DLANode*)parent)->DataType(node->InputConnNum(0)).m_dist, ddla->DataType(0).m_dist))) {
    return false;
  }
  return true;
}

void RemoveNOPRedistribs::Apply(Node *node) const
{
  Node *parent = node->Input(0);
  node->RedirectChildren(parent,node->InputConnNum(0));
  node->m_poss->DeleteChildAndCleanUp(node);
}

bool CombineRedistribs::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const Node *parent = node->Input(0);
  if (!parent)
    return false;
  NodeConnVecConstIter iter = parent->m_children.begin();
  for( ; iter != parent->m_children.end(); ++iter) {
    DLANode *output = (DLANode*)((*iter)->m_n);
    if (output != node
        && (output->GetNodeClass() == RedistNode::GetClass())
        && (output->InputConnNum(0) == node->InputConnNum(0))
        && DistTypeEqual(((RedistNode*)output)->m_info.m_dist, redist->m_info.m_dist)) 
      {
	return true;
      }
  }
  return false;
}

void CombineRedistribs::Apply(Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const Node *parent = node->Input(0);
  if (!parent)
    throw;
  NodeConnVecConstIter iter = parent->m_children.begin();
  for( ; iter != parent->m_children.end(); ++iter) {
    DLANode *output = (DLANode*)((*iter)->m_n);
    if (output != node
        && (output->GetNodeClass() == RedistNode::GetClass())
        && (output->InputConnNum(0) == node->InputConnNum(0))
        && DistTypeEqual(((RedistNode*)output)->m_info.m_dist, redist->m_info.m_dist)) 
      {
	output->RedirectChildren(node,0);
	output->m_poss->DeleteChildAndCleanUp(output);
	return;
      }
  }
  throw;
}

bool SplitRedistribs::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &src = redist->InputDataType(0).m_dist;
  const DistType *dest = &(redist->m_info.m_dist);

  if (src.m_numDims != dest->m_numDims)
    throw;
  else if (src.m_numDims <= m_dim)
    return false;
  else {
    if (src.m_dists[m_dim] != dest->m_dists[m_dim]) {
      //      cout << "src " << src.m_dists[m_dim] << endl;
      //      cout << "dest " << dest->m_dists[m_dim] << endl;
      DimVec destDims =  dest->m_dists[m_dim].DistEntryDims();
      for (Dim dim = 0; dim < dest->m_numDims; ++dim) {
	if (dim != m_dim) {
	  DistEntry srcDistEntry = src.m_dists[dim];
	  //	    unsigned int destDistEntry = dest->m_dists[dim];
	  DimVec srcDims = srcDistEntry.DistEntryDims();
	  DimVecIter iter = srcDims.begin();
	  for(; iter != srcDims.end(); ++iter) {
	    DimVecIter iter2 = destDims.begin();
	    for (; iter2 != destDims.end(); ++iter2) {
	      if (*iter == *iter2)
		return false;
	    }
	  }
	  /*
	    if (srcDistEntry != destDistEntry) {
	    if (srcDistEntry && !IsPrefix(srcDims,
	    DistType::DistEntryDims(destDistEntry)))
	    {
	    return true;
	    }
	    }
	  */
	}
      }
      return true;//false;
    }
    else {
      return false;
    }
  }
}

void SplitRedistribs::Apply(Node *node) const
{
  RedistNode *orig = (RedistNode*)node;
  DistType one = orig->InputDataType(0).m_dist;
  const DistType *two = &(orig->m_info.m_dist);

  one.m_dists[m_dim] = two->m_dists[m_dim];

/*  if (!one.IsSane()) {
    cout << "splitting " << m_dim << endl;
    cout << DistTypeToStr(one) << endl;
    cout << "from " << DistTypeToStr(orig->InputDataType(0).m_dist)
	 << " -> " << DistTypeToStr(orig->m_info.m_dist);
    throw;
    }*/



  RedistNode *newRedist = new RedistNode(one);
  newRedist->AddInput(orig->Input(0), orig->InputConnNum(0));
  node->m_poss->AddNode(newRedist);

  RedistNode *newRedist2 = new RedistNode(*two);
  newRedist2->AddInput(newRedist, 0);
  node->m_poss->AddNode(newRedist2);

  node->RedirectChildren(newRedist2, 0);
  node->m_poss->DeleteChildAndCleanUp(node);

}

bool SingleIndexAllToAll::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).m_dist;

  if (srcType.m_numDims != redist->m_info.m_dist.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_dim)
    return false;

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = redist->m_info.m_dist.m_dists[m_dim];
  if (srcEntry == destEntry)
    return false;

  for(Dim dim = 0; dim < srcType.m_numDims; ++dim) {
    if (dim != m_dim) {
      if (srcType.m_dists[dim] != redist->m_info.m_dist.m_dists[dim])
	return false;
    }
  }

  if (srcEntry.IsStar() || destEntry.IsStar())
    return false;
  
  DimVec srcDims = srcEntry.DistEntryDims();
  DimVec destDims = destEntry.DistEntryDims();

  if (IsPrefix(srcDims,destDims) || IsPrefix(destDims,srcDims))
    return false;

  if (srcDims.size() == destDims.size()) {
    DimSet srcSet;
    srcSet.insert(srcDims.begin(), srcDims.end());
    DimSet destSet;
    destSet.insert(destDims.begin(), destDims.end());
    if (includes(srcSet.begin(), srcSet.end(),
		 destSet.begin(), destSet.end()))
      return false;
  }

  return true;
}

void SingleIndexAllToAll::Apply(Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).m_dist;

  DistType type1 = srcType;
  DistEntry entry1 = type1.m_dists[m_dim];
  DimVec vec1 = entry1.DistEntryDims();
  DimSet set1;
  set1.insert(vec1.begin(), vec1.end());
  
  DimVec destVec = redist->m_info.m_dist.m_dists[m_dim].DistEntryDims();
  DimVecIter iter = destVec.begin();
  for(; iter != destVec.end(); ++iter) {
    Dim dim = *iter;
    if (set1.find(dim) == set1.end()) {
      vec1.push_back(dim);
    }
  }
  type1.m_dists[m_dim].DimsToDistEntry(vec1);

  RedistNode *redist1 = new RedistNode(type1);
  redist1->AddInput(redist->Input(0), redist->InputConnNum(0));
  node->m_poss->AddNode(redist1);

  DistType type2 = redist->m_info.m_dist;
  DistEntry entry2 = type2.m_dists[m_dim];
  DimVec vec2 = entry2.DistEntryDims();
  DimSet set2;
  set2.insert(vec2.begin(), vec2.end());
  
  DimVec srcVec = srcType.m_dists[m_dim].DistEntryDims();
  iter = srcVec.begin();
  for(; iter != srcVec.end(); ++iter) {
    Dim dim = *iter;
    if (set2.find(dim) == set2.end()) {
      vec2.push_back(dim);
    }
  }
  type2.m_dists[m_dim].DimsToDistEntry(vec2);
  
  RedistNode *redist2 = new RedistNode(type2);
  redist2->AddInput(redist1, 0);
  node->m_poss->AddNode(redist2);

  RedistNode *redist3 = new RedistNode(redist->m_info.m_dist);
  redist3->AddInput(redist2, 0);
  node->m_poss->AddNode(redist3);

  node->RedirectChildren(redist3, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}


bool DoubleIndexAllToAll::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).m_dist;

  if (srcType.m_numDims != redist->m_info.m_dist.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_dim)
    return false;

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = redist->m_info.m_dist.m_dists[m_dim];
  if (srcEntry == destEntry)
    return false;

  for(Dim dim = m_dim+1; dim < srcType.m_numDims; ++dim) {
    DistEntry srcEntry2 = srcType.m_dists[dim];
    DistEntry destEntry2 = redist->m_info.m_dist.m_dists[dim];
    if (srcEntry2 != destEntry2) {
      DimVec srcVec = srcEntry.DistEntryDims();
      DimVec destVec = destEntry.DistEntryDims();
      DimVec srcVec2 = srcEntry2.DistEntryDims();
      DimVec destVec2 = destEntry2.DistEntryDims();

      //[(x|_|y),(z|_|w)] <- [(x|_|u|_|w),(z|_|v|_|y)]

      DimVec suff1, suff2;
      DimVec pref11, pref12, pref21, pref22;
      // suff1 = w
      // pref11 = x|_|u
      // pref12 = z
      GetCommonSuffix(srcVec, destVec2, suff1, pref11, pref12);
      // suff2 = y
      // pref21 = x
      // pref22 = z|_|v
      GetCommonSuffix(destVec, srcVec2, suff2, pref21, pref22);
      if (!(suff1.empty()&&suff2.empty()))  {
	if ((pref11 == pref21 || pref21.empty()) 
	    && (pref12 == pref22 || pref12.empty())) {
	  DistType intType = srcType;

	  DimVec newDims1;
	  DimVec newDims2;
	    
		  
	  // have to choose where u and v go

	  if (!suff1.empty() && !suff2.empty()) {
	    //[(x|_|u|_|y),(z|_|v|_|w)] <- [(x|_|u|_|w),(z|_|v|_|y)]
	    //[(x|_|y),(z|_|w)] <- [(x|_|u|_|y),(z|_|v|_|w)]
	    newDims1 = pref11;
	    newDims1.insert(newDims1.end(),
			    suff2.begin(), suff2.end());
	    newDims2 = pref22;
	    newDims2.insert(newDims2.end(),
			    suff1.begin(), suff1.end());		
	  }
	  else {
	    if (suff1.empty())  { // w is empty, y isn't
	      //[(x|_|y),(z|_|v|_|u)] <- [(x|_|u),(z|_|v|_|y)]
	      //[(x|_|y),(z)] <- [(x|_|y),(z|_|v|_|u)]
	      newDims1 = GetCommonPrefix(srcVec, destVec); // x
	      
	      newDims2 = pref22; // z |_| v
	      newDims2.insert(newDims2.end(), 
			      srcVec.begin()+newDims1.size(), srcVec.end()); // add u
	      
	      newDims1.insert(newDims1.end(), suff2.begin(), suff2.end()); // add y
	    }
	    else { // y is empty, w isn't
	      //[(x|_|u|_|v),(z|_|w)] <- [(x|_|u|_|w),(z|_|v)]
	      //[(x),(z|_|w)] <- [(x|_|u|_|v),(z|_|w)]
	      newDims2 = GetCommonPrefix(srcVec2, destVec2); // z
	      
	      newDims1 = pref11; // x |_| u
	      newDims1.insert(newDims1.end(), 
			      srcVec2.begin()+newDims2.size(), srcVec2.end()); // add v
	      
	      newDims2.insert(newDims2.end(), suff1.begin(), suff1.end()); // add w
	    }
	  }	    
	  intType.m_dists[m_dim].DimsToDistEntry(newDims1);
	  intType.m_dists[dim].DimsToDistEntry(newDims2);
	  if (DistTypeNotEqual(intType, redist->m_info.m_dist))
	    return true;
	}
      }
    }
  }
  return false;
}

void DoubleIndexAllToAll::Apply(Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).m_dist;

  if (srcType.m_numDims != redist->m_info.m_dist.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_dim)
    throw;

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = redist->m_info.m_dist.m_dists[m_dim];
  if (srcEntry == destEntry)
    throw;

  for(Dim dim = m_dim+1; dim < srcType.m_numDims; ++dim) {
    DistEntry srcEntry2 = srcType.m_dists[dim];
    DistEntry destEntry2 = redist->m_info.m_dist.m_dists[dim];
    if (srcEntry2 != destEntry2) {
      DimVec srcVec = srcEntry.DistEntryDims();
      DimVec destVec = destEntry.DistEntryDims();
      DimVec srcVec2 = srcEntry2.DistEntryDims();
      DimVec destVec2 = destEntry2.DistEntryDims();

      //[(x|_|y),(z|_|w)] <- [(x|_|u|_|w),(z|_|v|_|y)]

      DimVec suff1, suff2;
      DimVec pref11, pref12, pref21, pref22;
      // suff1 = w
      // pref11 = x|_|u
      // pref12 = z
      GetCommonSuffix(srcVec, destVec2, suff1, pref11, pref12);
      // suff2 = y
      // pref21 = x
      // pref22 = z|_|v
      GetCommonSuffix(destVec, srcVec2, suff2, pref21, pref22);

      if (!(suff1.empty()&&suff2.empty()))  {
	if ((pref11 == pref21 || pref21.empty()) 
	    && (pref12 == pref22 || pref12.empty())) {
	  DistType intType = srcType;

	  DimVec newDims1;
	  DimVec newDims2;
	    
		  
	  // have to choose where u and v go

	  if (!suff1.empty() && !suff2.empty()) {
	    //[(x|_|u|_|y),(z|_|v|_|w)] <- [(x|_|u|_|w),(z|_|v|_|y)]
	    //[(x|_|y),(z|_|w)] <- [(x|_|u|_|y),(z|_|v|_|w)]
	    newDims1 = pref11;
	    newDims1.insert(newDims1.end(),
			    suff2.begin(), suff2.end());
	    newDims2 = pref22;
	    newDims2.insert(newDims2.end(),
			    suff1.begin(), suff2.end());		
	  }
	  else {

	    if (suff1.empty())  { // w is empty, y isn't
	      //[(x|_|y),(z|_|v|_|u)] <- [(x|_|u),(z|_|v|_|y)]
	      //[(x|_|y),(z)] <- [(x|_|y),(z|_|v|_|u)]
	      newDims1 = GetCommonPrefix(srcVec, destVec); // x
	      
	      newDims2 = pref22; // z |_| v
	      newDims2.insert(newDims2.end(), 
			      srcVec.begin()+newDims1.size(), srcVec.end()); // add u
	      
	      newDims1.insert(newDims1.end(), suff2.begin(), suff2.end()); // add y
	    }
	    else { // y is empty, w isn't
	      //[(x|_|u|_|v),(z|_|w)] <- [(x|_|u|_|w),(z|_|v)]
	      //[(x),(z|_|w)] <- [(x|_|u|_|v),(z|_|w)]
	      newDims2 = GetCommonPrefix(srcVec2, destVec2); // z
	      
	      newDims1 = pref11; // x |_| u
	      newDims1.insert(newDims1.end(), 
			      srcVec2.begin()+newDims2.size(), srcVec2.end()); // add v
	      
	      newDims2.insert(newDims2.end(), suff1.begin(), suff1.end()); // add w
	    }
	  }	    


	  intType.m_dists[m_dim].DimsToDistEntry(newDims1);
	  intType.m_dists[dim].DimsToDistEntry(newDims2);
	  if (DistTypeNotEqual(intType, redist->m_info.m_dist)) {
	      RedistNode *newRedist = new RedistNode(intType);
	      newRedist->AddInput(redist->Input(0), redist->InputConnNum(0));
	      node->m_poss->AddNode(newRedist);
	    
	      redist->ChangeInput2Way(redist->Input(0), redist->InputConnNum(0),
				      newRedist, 0);    
	  
	      return;	    
	  }
	}
      }
    }
  }
  throw;
}

bool SplitAllGathers::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).m_dist;
  const DistType &destType = redist->m_info.m_dist;

  if (srcType.m_numDims != redist->m_info.m_dist.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_dim)
    return false;

  for (Dim dim = 0; dim < srcType.m_numDims; ++dim) {
    if (dim != m_dim) {
      if (srcType.m_dists[dim] != destType.m_dists[dim])
	return false;
    }
  }

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = destType.m_dists[m_dim];
  if (srcEntry == destEntry)
    return false;

  DimVec srcDims = srcEntry.DistEntryDims();

  if (destEntry.IsStar() && srcDims.size() > 1)
    return true;

  DimVec destDims = destEntry.DistEntryDims();

  if (srcDims.size() < destDims.size())
    return false;

  if ((srcDims.size() - destDims.size()) > 1 &&  IsPrefix(destDims, srcDims))
    return true;
  else
    return false;
}

void SplitAllGathers::Apply(Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).m_dist;
  const DistType &destType = redist->m_info.m_dist;

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = destType.m_dists[m_dim];

  DimVec destDims = destEntry.DistEntryDims();
  DimVec srcDims = srcEntry.DistEntryDims();

  srcDims.pop_back();

  DistType intType = srcType;
  intType.m_dists[m_dim].DimsToDistEntry(srcDims);

  RedistNode *redist1 = new RedistNode(intType);
  redist1->AddInput(redist->Input(0), redist->InputConnNum(0));
  node->m_poss->AddNode(redist1);

  RedistNode *redist2 = new RedistNode(destType);
  redist2->AddInput(redist1, 0);
  node->m_poss->AddNode(redist2);

  node->RedirectChildren(redist2, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

bool CombineDisappearingModes::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).m_dist;
  const DistType &destType = redist->m_info.m_dist;

  if (srcType.m_numDims != redist->m_info.m_dist.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_srcDim || srcType.m_numDims <= m_destDim)
    return false;

  DimSet usedDims = destType.UsedGridDims();

  DimVec dest = srcType.m_dists[m_destDim].DistEntryDims();
  DimVec src = srcType.m_dists[m_srcDim].DistEntryDims();
  
  if (dest.empty() || src.empty())
    return false;
  if (usedDims.find(dest.back()) == usedDims.end()) {
    return usedDims.find(src.back()) == usedDims.end();
  }
  return false;
}

void CombineDisappearingModes::Apply(Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).m_dist;
  const DistType &destType = redist->m_info.m_dist;

  if (srcType.m_numDims != redist->m_info.m_dist.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_srcDim || destType.m_numDims <= m_destDim)
    throw;

  DistType intType = srcType;
  DimVec dest = intType.m_dists[m_destDim].DistEntryDims();
  DimVec src = intType.m_dists[m_srcDim].DistEntryDims();
  dest.push_back(src.back());
  src.pop_back();
  intType.m_dists[m_destDim].DimsToDistEntry(dest);
  intType.m_dists[m_srcDim].DimsToDistEntry(src);

  RedistNode *newRedist = new RedistNode(intType);
  newRedist->AddInput(redist->Input(0), redist->InputConnNum(0));
  redist->m_poss->AddNode(newRedist);
  
  redist->ChangeInput2Way(redist->Input(0), redist->InputConnNum(0),
			  newRedist, 0);    
}

bool PermuteDistribution::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;

  const DistType &srcType = redist->InputDataType(0).m_dist;
  const DistType &destType = redist->m_info.m_dist;

  if (srcType.m_numDims != redist->m_info.m_dist.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_srcDim || srcType.m_numDims <= m_destDim)
    return false;


  const DistEntry &srcEntry = srcType.m_dists[m_srcDim];
  const DistEntry &destEntry = destType.m_dists[m_destDim];

  if (srcEntry.IsStar() || destEntry.IsStar())
    return false;

  DimSet usedSrcDims = srcEntry.DistEntryDimSet();
  DimSet usedDestDims = destEntry.DistEntryDimSet();

  DimSet intersection;
  std::set_intersection(usedSrcDims.begin(), usedSrcDims.end(),
			usedDestDims.begin(), usedDestDims.end(),
			std::inserter(intersection, intersection.begin()));
  
  DimVec vec = srcEntry.DistEntryDims();
  DimVecIter iter = vec.begin();
  do {
    if (intersection.find(*iter) != intersection.end()) {
      vec.erase(iter);
      iter = vec.begin();
    }
    else
      ++iter;
  } while (iter != vec.end());

  DimVec destVec = destEntry.DistEntryDims();
  iter = destVec.begin();
  for(; iter != destVec.end(); ++iter) {
    Dim dim = *iter;
    if (intersection.find(dim) != intersection.end())
      vec.push_back(dim);
  }

  DistType type = srcType;
  type.m_dists[m_srcDim].DimsToDistEntry(vec);

  return DistTypeNotEqual(type, srcType);
}

void PermuteDistribution::Apply(Node *node) const
{  
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  RedistNode *redist = (RedistNode*)node;

  const DistType &srcType = redist->InputDataType(0).m_dist;
  const DistType &destType = redist->m_info.m_dist;

  if (srcType.m_numDims != redist->m_info.m_dist.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_srcDim || srcType.m_numDims <= m_destDim)
    throw;


  const DistEntry &srcEntry = srcType.m_dists[m_srcDim];
  const DistEntry &destEntry = destType.m_dists[m_destDim];

  DimSet usedSrcDims = srcEntry.DistEntryDimSet();
  DimSet usedDestDims = destEntry.DistEntryDimSet();

  DimSet intersection;
  std::set_intersection(usedSrcDims.begin(), usedSrcDims.end(),
			usedDestDims.begin(), usedDestDims.end(),
			std::inserter(intersection, intersection.begin()));
  
  DimVec vec = srcEntry.DistEntryDims();
  DimVecIter iter = vec.begin();
  do {
    if (intersection.find(*iter) != intersection.end()) {
      vec.erase(iter);
      iter = vec.begin();
    }
    else
      ++iter;
  } while (iter != vec.end());

  DimVec destVec = destEntry.DistEntryDims();
  iter = destVec.begin();
  for(; iter != destVec.end(); ++iter) {
    if (intersection.find(*iter) != intersection.end())
      vec.push_back(*iter);
  }

  DistType type = srcType;
  type.m_dists[m_srcDim].DimsToDistEntry(vec);

  RedistNode *newRedist = new RedistNode(type);
  newRedist->AddInput(redist->Input(0), redist->InputConnNum(0));
  redist->m_poss->AddNode(newRedist);
  
  redist->ChangeInput2Way(redist->Input(0), redist->InputConnNum(0),
			  newRedist, 0);    
}


#endif



