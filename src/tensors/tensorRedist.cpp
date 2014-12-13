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

bool GetAllToAllPattern(const DistType &srcType, 
			const DistType &destType,
			DistEntryVec *gridModesInvolved);

bool GetAllGatherPattern(const DistType &srcType,
			 const DistType &destType,
			 DistEntryVec *gridModes);

bool GetLocalRedistPattern(const DistType &srcType,
			   const DistType &destType,
			   DimVec *indices,
			   DistEntryVec *gridModes);
			

RedistNode::RedistNode() 
{
  m_info.SetToDefault(0);
  m_lsizes = NULL;
}

RedistNode::RedistNode(const DistType &destType, const string &align, 
	     const DimVec &alignModes, const DimVec &alignModesSrc)
{
  m_info.SetDistAndClearPerm(destType);
  m_lsizes = NULL;
  m_align = align;
  m_alignModes = alignModes;
  m_alignModesSrc = alignModesSrc;
}

RedistNode::RedistNode(const DistType &destType, const Permutation &perm, const string &align, 
	     const DimVec &alignModes, const DimVec &alignModesSrc)
  : RedistNode(destType, align, alignModes, alignModesSrc)
{
  m_info.SetPerm(perm);
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
  m_align = origNode->m_align;
  m_alignModes = origNode->m_alignModes;
  m_alignModesSrc = origNode->m_alignModesSrc;
}

NodeType RedistNode::GetType() const
{
  if (m_name.length()) 
    return m_name;
  else {
    return (string)"RedistNode to " +  m_info.GetDist().QuickStr() + m_info.GetPerm().Str();
  }
}

void RedistNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLANode::Prop();

    if (m_align.empty())
      throw;

    if (m_inputs.size() != 1) {
      cout << "m_inputs.size() != 1\n";
      throw;
    }

    if (!m_children.size())
      throw;
  
    if (!m_name.length())
      m_name = (string)"RedistNode to " +  m_info.GetDist().QuickStr() + m_info.GetPerm().Str();
    DLANode *parent = (DLANode*)Input(0);
    parent->Prop();

    //if the input has a permutation, then some of the below might need to be 
    // redone
    if (parent->DataType(InputConnNum(0)).HasPerm())
      throw;


    const DistType &m_srcType = parent->DataType(InputConnNum(0)).GetDist();
    const Dim numDims = m_info.GetDist().m_numDims;

    if (m_info.GetDist() == m_srcType)
      throw;

    if (m_info.HasPerm() && m_info.GetPerm().Size() != numDims)
      throw;

    if (numDims != InputNumDims(0)) {
      cout << "numDims " << numDims << " " << InputNumDims(0) << endl;
      cout << "input " << Input(0)->GetNodeClass() << endl;
      cout << GetInputNameStr(0) << endl;
      cout << m_info.GetDist().PrettyStr() << endl;
      throw;
    }
    if (m_srcType.m_numDims != numDims)
      throw;

    if (m_srcType.m_numDims != numDims) {
      cout << m_info.GetDist().str() << " <- " << m_srcType.str() << endl;
      throw;
    }

    if (!m_info.GetDist().IsSane()) {
      cout << m_info.GetDist().str() << endl;
      m_poss->PrintTransVec();
      cout << m_srcType.PrettyStr() << endl;
      cout <<m_info.GetDist().PrettyStr() << endl;
      throw;
    }

    m_cost = 0;

    DistEntryVec gridModesInvolved;
    if (GetAllToAllPattern(m_srcType, m_info.GetDist(),
			   &gridModesInvolved)) {
      unsigned int numProcs = 1;
      DistEntryVecIter gridModeIter = gridModesInvolved.begin();
      for(; gridModeIter != gridModesInvolved.end(); ++gridModeIter) {
	DistEntry entry = *gridModeIter;
	DimVec gridModes = entry.DistEntryDims();
	DimVecIter iter = gridModes.begin();
	for(; iter != gridModes.end(); ++iter) {
	  numProcs *= GridLens[*iter];
	}
      }
      const unsigned int totNumIters = m_lsizes[0].NumSizes();
      for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
	Cost temp = 1;
	for (Dim dim = 0; dim < numDims; ++dim) {
	  temp *= m_lsizes[dim][iteration];
	}
	m_cost += AllToAll(temp, numProcs);
	m_cost += (PSIR+PSIW)*(2*temp);
      }
      return;
    }

    if (GetAllGatherPattern(m_srcType, m_info.GetDist(),
			    &gridModesInvolved)) {
      unsigned int numProcs = 1;
      DistEntryVecIter gridModeIter = gridModesInvolved.begin();
      for(; gridModeIter != gridModesInvolved.end(); ++gridModeIter) {
	DistEntry entry = *gridModeIter;
	DimVec gridModes = entry.DistEntryDims();
	DimVecIter iter = gridModes.begin();
	for(; iter != gridModes.end(); ++iter) {
	  numProcs *= GridLens[*iter];
	}
      }
      /*
      cout << "***\n";
      cout << "For " << InputDataType(0).m_dist.PrettyStr() << " -> "
	   << m_info.m_dist.PrettyStr() << endl;
      */
      const unsigned int totNumIters = m_lsizes[0].NumSizes();
      for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
	Cost temp = 1;
	for (Dim dim = 0; dim < numDims; ++dim) {
	  temp *= m_lsizes[dim][iteration];
	}
	//	cout << "AllGather( " << std::scientific << temp * numProcs << ", " << numProcs << " )\n";
	//	cout << "\t" << temp << " data\n";
	m_cost += AllGather(temp, numProcs);
	m_cost += (PSIR+PSIW)*(temp + temp / numProcs);
	//	cout << "cost " << m_cost << endl;
      }
      //      cout << "***\n";
      return;
    }

    DimVec indices;
    if (GetLocalRedistPattern(m_srcType, m_info.GetDist(),
			      &indices,
			      &gridModesInvolved)) {
      //Local mem copy
      m_cost = (PSIR+PSIW)*(TotalNumberOfLocalElements(0));
      return;
    }

    DimSet diffs;
  
    for (Dim dim = 0; dim < numDims; ++dim) {
      if (m_srcType.m_dists[dim] != m_info.GetDist().m_dists[dim]) {
	diffs.insert(dim);
      }
    }
  
    if (diffs.empty()) {
      throw;
    }
    else if (diffs.size() == 1) {
      Dim dim = *(diffs.begin());
      DimVec src = m_srcType.m_dists[dim].DistEntryDims();
      DimVec dest = m_info.GetDist().m_dists[dim].DistEntryDims();

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
	m_cost += AllToAll(temp, numProcs);
	m_cost += (PSIR+PSIW)*(2*temp);
      }
    }
    else
      m_cost = 0;
  
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
  /*
  if (iter1 == dims1.end())
    throw;
  if (iter2 == dims2.end())
    throw;
  */
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
  if (iter2 == dims2.end()) {
    DistEntry entry1, entry2;
    entry1.DimsToDistEntry(dims1);
    cout << entry1.PrettyStr() << endl;
    entry2.DimsToDistEntry(dims2);
    cout << entry2.PrettyStr() << endl;
    throw;
  }
  while (iter2 != dims2.end()) {
    suff.push_back(*iter2);
    ++iter2;
  }
}


Phase RedistNode::MaxPhase() const
{
  DLANode *parent = (DLANode*)Input(0);
  const DistType &m_srcType = parent->DataType(InputConnNum(0)).GetDist();
  const Dim numDims = m_info.GetDist().m_numDims;
  

  //Check for multi-mode AllToAll

  if (IsPrimitive()) 
    {
      return NUMPHASES;
    }
  
  DimSet diffs;
  
  for (Dim dim = 0; dim < numDims; ++dim) {
    if (m_srcType.m_dists[dim] != m_info.GetDist().m_dists[dim]) {
      diffs.insert(dim);
    }
  }
  
  if (diffs.empty()) {
    throw;
  }
  else if (diffs.size() == 1) {
    Dim dim = *(diffs.begin());
    DimVec src = m_srcType.m_dists[dim].DistEntryDims();
    DimVec dest = m_info.GetDist().m_dists[dim].DistEntryDims();
    if (m_srcType.m_dists[dim].DistEntryDimSet() != m_info.GetDist().m_dists[dim].DistEntryDimSet()) {
      return ROTENSORPHASE;
    }
    else
      return NUMPHASES;
  }
  else
    return ROTENSORPHASE;
}

bool RedistNode::IsPrimitive() const
{
  const DLANode *parent = (DLANode*)Input(0);
  const DistType &m_srcType = parent->DataType(InputConnNum(0)).GetDist();
  const Dim numDims = m_info.GetDist().m_numDims;
  

  //Check for multi-mode AllToAll

  if (GetAllToAllPattern(m_srcType, m_info.GetDist(),
			 NULL)) 
    {
      return true;
    }

  if (GetAllGatherPattern(m_srcType, m_info.GetDist(),
			  NULL))
    {
      return true;
    }

  if (GetLocalRedistPattern(m_srcType, m_info.GetDist(),
			    NULL, NULL))
    {
      return true;
    }
  
  DimSet diffs;
  
  for (Dim dim = 0; dim < numDims; ++dim) {
    if (m_srcType.m_dists[dim] != m_info.GetDist().m_dists[dim]) {
      diffs.insert(dim);
    }
  }
  
  if (diffs.empty()) {
    throw;
  }
  else if (diffs.size() == 1) {
    Dim dim = *(diffs.begin());
    DimVec src = m_srcType.m_dists[dim].DistEntryDims();
    DimVec dest = m_info.GetDist().m_dists[dim].DistEntryDims();
    if (m_srcType.m_dists[dim].DistEntryDimSet() != m_info.GetDist().m_dists[dim].DistEntryDimSet()) {
      return false;
    }
    return true;
  }
  else
    return false;
}

const Dim RedistNode::NumDims(ConnNum num) const
{
  if (num > 0)
    throw;
  return InputNumDims(0);
}

const Sizes* RedistNode::Len(ConnNum num, Dim dim) const
{
  if (num > 0)
    throw;
  if (InputDataType(0).HasPerm())
    throw;
  if (m_info.HasPerm())
    return InputLen(0,m_info.GetPerm().MapFinishToStart(dim));
  else
    return InputLen(0,dim);
}

const Sizes* RedistNode::LocalLen(ConnNum num, Dim dim) const
{
  if (num > 0)
    throw;
  if (InputDataType(0).HasPerm())
    throw;
  if (!m_isArray) {
    return m_lsizes;
  }
  else {
    if (m_info.HasPerm())
      return m_lsizes+m_info.GetPerm().MapFinishToStart(dim);
    else
      return m_lsizes+dim;
  }
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

  //  cout << "For " << InputDataType(0).m_dist.PrettyStr() << " -> "
  //       << m_info.m_dist.PrettyStr() << endl;

  DLANode *in = (DLANode*)Input(0);
  ConnNum num = InputConnNum(0);
  Dim numDims = in->NumDims(num);
  if (numDims) {
    m_isArray = true;
    m_lsizes = new Sizes[numDims];
    for (Dim dim = 0; dim < numDims; ++dim) {
      GetLocalSizes(m_info.GetDist(), dim, in->Len(num,dim), m_lsizes+dim); 
      //      cout << "dim " << dim << ": ";
      //      in->Len(num,dim)->Print();
      //      cout << "*to*\n";
      //      (m_lsizes+dim)->Print();	
    }
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

Name RedistNode::GetName(ConnNum num) const
{
  if (num > 0)
    throw;
  Name name = GetInputName(0);
  name.m_permutation = m_info.GetPerm();
  name.m_type = m_info.GetDist();
  return name;
}

const DataTypeInfo& RedistNode::DataType(ConnNum num) const 
{
  return m_info;
}


void RedistNode::PrintCode(IndStream &out)
{  
  //Reflect in AddVars
  const DistType &m_srcType = InputDataType(0).GetDist();
  const Dim numDims = m_info.GetDist().m_numDims;

  string inName = GetInputName(0).str();
  string outName = GetName(0).str(); 

  out.Indent();
  *out << "   // " << GetName(0).PrettyStr() 
       << " <- " << GetInputName(0).PrettyStr();
  if (InputDataType(0).HasPerm())
    *out << " with permutation " << InputDataType(0).GetPerm().Str();
  *out << endl;


  if (m_align.empty()) 
    throw;
  if (m_alignModes.size() != m_alignModesSrc.size())
    throw;


  bool align = false;
  if (m_alignModes.empty())
    align = false; // to superscede any of the following
  else if (m_align != GetInputNameStr(0)) 
    align = true;
  else {
    //aligning to input, so make sure the alignment isn't the identity
    if (m_alignModes != m_alignModesSrc)
      align = true;
  }

  if (align) {
    out.Indent();
    *out << outName << ".AlignModesWith( " 
	 << ModeArrayVarName(m_alignModes) << ", "
	 << m_align << ", " 
	 << ModeArrayVarName(m_alignModesSrc) << " );\n";
  }


  out.Indent();

  if (CurrPhase <= ROTENSORPHASE)
    return;
  
  if (m_srcType.m_numDims != numDims)
    throw;

  DistEntryVec gridModesInvolved;
  if (GetAllToAllPattern(m_srcType, m_info.GetDist(),
			 &gridModesInvolved)) {
    *out << outName << ".AllToAllRedistFrom( " 
	 << inName << ", "
	 << DistEntryVecVarName(gridModesInvolved) << " );\n";
    return;
  }

  DimVec indices;
  if (GetAllGatherPattern(m_srcType, m_info.GetDist(),
			  &gridModesInvolved)) {
    *out << outName << ".AllGatherRedistFrom( "
	 << inName << ", "
	 << DistEntryVecVarName(gridModesInvolved) << " );\n";
    return;
  }

  if (GetLocalRedistPattern(m_srcType, m_info.GetDist(),
			  &indices,
			  &gridModesInvolved)) {
    *out << outName << ".LocalRedistFrom( "
	 << inName << ", "
      	 << ModeArrayVarName(indices) << ", "
	 << DistEntryVecVarName(gridModesInvolved) << " );\n";
    return;
  }

  DimSet diffs;

  for (Dim dim = 0; dim < numDims; ++dim) {
    if (m_srcType.m_dists[dim] != m_info.GetDist().m_dists[dim]) {
      diffs.insert(dim);
    }
  }

  if (diffs.empty()) {
    throw;
  }
  else if (diffs.size() == 1) {
    Dim dim = *(diffs.begin());
    DimVec src = m_srcType.m_dists[dim].DistEntryDims();
    DimVec dest = m_info.GetDist().m_dists[dim].DistEntryDims();
    if (m_srcType.m_dists[dim].DistEntryDimSet() != m_info.GetDist().m_dists[dim].DistEntryDimSet())
      throw;       
    if (m_srcType.m_dists[dim].DistEntryDimSet() != m_info.GetDist().m_dists[dim].DistEntryDimSet())
      throw;
    *out << outName << ".PermutationRedistFrom( "
	 << inName << ", " << dim << ", "
	 << ModeArrayVarName(src) << " );\n";
    return;
  }
  else 
    throw;
}

void RedistNode::AddVariables(VarSet &set) const
{
  //Reflect in PrintCode
  DLANode::AddVariables(set);

  if (CurrPhase <= ROTENSORPHASE)
    return;

  {
    Var var(ModeArrayVarType, m_alignModes);
    set.insert(var);
  }
  {
    Var var(ModeArrayVarType, m_alignModesSrc);
    set.insert(var);
  }


  
  const DistType &m_srcType = InputDataType(0).GetDist();
  const Dim numDims = m_info.GetDist().m_numDims;

  string inName = GetInputName(0).str();
  string outName = GetName(0).str(); 

  if (m_srcType.m_numDims != numDims)
    throw;

  DistEntryVec gridModesInvolved;
  if (GetAllToAllPattern(m_srcType, m_info.GetDist(),
			 &gridModesInvolved)) {
    { 
      DistEntryVecIter iter = gridModesInvolved.begin();
      for(; iter != gridModesInvolved.end(); ++iter) {
	Var var(ModeArrayVarType, (*iter).DistEntryDims());
	set.insert(var);
      }
      Var var(gridModesInvolved);
      set.insert(var);
    }

    return;
  }

  if (GetAllGatherPattern(m_srcType, m_info.GetDist(),
			  &gridModesInvolved)) {
    { 
      DistEntryVecIter iter = gridModesInvolved.begin();
      for(; iter != gridModesInvolved.end(); ++iter) {
	Var var(ModeArrayVarType, (*iter).DistEntryDims());
	set.insert(var);
      }
      Var var(gridModesInvolved);
      set.insert(var);
    }
    return;
  }

  DimVec indices;
  if (GetLocalRedistPattern(m_srcType, m_info.GetDist(),
			  &indices,
			  &gridModesInvolved)) {
    {
      Var var(ModeArrayVarType, indices);
      set.insert(var);
    }
    { 
      DistEntryVecIter iter = gridModesInvolved.begin();
      for(; iter != gridModesInvolved.end(); ++iter) {
	Var var(ModeArrayVarType, (*iter).DistEntryDims());
	set.insert(var);
      }
      Var var(gridModesInvolved);
      set.insert(var);
    }
    return;
  }

  DimSet diffs;
    

  for (Dim dim = 0; dim < numDims; ++dim) {
    if (m_srcType.m_dists[dim] != m_info.GetDist().m_dists[dim]) {
      diffs.insert(dim);
    }
  }

  if (diffs.empty()) {
    cout << m_info.GetDist().PrettyStr() << " <- " << m_srcType.PrettyStr() << endl;
    cout << inName << endl;
    throw;
  }
  else if (diffs.size() == 1) {
    Dim dim = *(diffs.begin());
    DimVec src = m_srcType.m_dists[dim].DistEntryDims();
    DimVec dest = m_info.GetDist().m_dists[dim].DistEntryDims();
    if (m_srcType.m_dists[dim].DistEntryDimSet() != m_info.GetDist().m_dists[dim].DistEntryDimSet()) {
      cout << m_srcType.PrettyStr() << " -> " << m_info.GetDist().PrettyStr() << endl;
      
      GetAllToAllPattern(m_srcType, m_info.GetDist(),
			 &gridModesInvolved);
      throw;
    }
    Var var(ModeArrayVarType, src);
    set.insert(var);
  }
  else
    throw;
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
    ConnNum num = InputConnNum(0);

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
  const DistType* type = &(redistNode->m_info.GetDist());
  if (node->m_children.size() == 0)
    throw;
  while(redistNode->Input(0) 
        && (redistNode->Input(0)->GetNodeClass() == RedistNode::GetClass()))
    {
      redistNode = (RedistNode*)redistNode->Input(0);
      if (DistTypeEqual(redistNode->m_info.GetDist(),*type))
	return true;
      for(ConnNum i = 0; i < redistNode->m_children.size(); ++i) {
	Node *tmp = redistNode->Child(i);
	if (tmp != node && tmp->GetNodeClass() == RedistNode::GetClass()) {
	  RedistNode *redist = ((RedistNode*)tmp);
	  if (DistTypeEqual(redist->m_info.GetDist(), *type))
	    return redist->m_info.GetPerm() == redistNode->m_info.GetPerm();
	}
      }
    }
  if (redistNode->Input(0)
      && DistTypeEqual(((DLANode*)(redistNode->Input(0)))->DataType(redistNode->InputConnNum(0)).GetDist(),*type))
    return true;
  return false;
}

void RemoveWastedRedist::Apply(Node *node) const
{
  RedistNode *redistNode = (RedistNode*)node;
  const DistType* type = &(redistNode->m_info.GetDist());
  while(redistNode->Input(0) 
        && (redistNode->Input(0)->GetNodeClass() == RedistNode::GetClass()))
    {
      redistNode = (RedistNode*)redistNode->Input(0);
      if (DistTypeEqual(redistNode->m_info.GetDist(), *type)) {
	node->RedirectChildren(redistNode, 0);
	node->m_poss->DeleteChildAndCleanUp(node);
	return;
      }
      for(ConnNum i = 0; i < redistNode->m_children.size(); ++i) {
	Node *tmp = redistNode->Child(i);
	if (tmp != node && tmp->GetNodeClass() == RedistNode::GetClass()) {
	  RedistNode *redist = ((RedistNode*)tmp);
	  if (DistTypeEqual(redist->m_info.GetDist(),*type)) {
	    if (redist->m_info.GetPerm() == redistNode->m_info.GetPerm()) {
	      node->RedirectChildren(tmp, 0);
	      node->m_poss->DeleteChildAndCleanUp(node);
	      return;
	    }
	  }
	}
      }
    }
  if (redistNode->Input(0)
      && DistTypeEqual(((DLANode*)redistNode->Input(0))->DataType(redistNode->InputConnNum(0)).GetDist(),*type))
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
      || (DistTypeNotEqual(((DLANode*)parent)->DataType(node->InputConnNum(0)).GetDist(), ddla->DataType(0).GetDist()))) {
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
        && DistTypeEqual(((RedistNode*)output)->m_info.GetDist(), redist->m_info.GetDist())
	&& ((RedistNode*)output)->m_info.GetPerm() == redist->m_info.GetPerm())
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
        && DistTypeEqual(((RedistNode*)output)->m_info.GetDist(), redist->m_info.GetDist())
	&& ((RedistNode*)output)->m_info.GetPerm() == redist->m_info.GetPerm())
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
  const DistType &src = redist->InputDataType(0).GetDist();
  const DistType *dest = &(redist->m_info.GetDist());

  if (src.m_numDims != dest->m_numDims)
    throw;
  else if (src.m_numDims <= m_dim)
    return false;
  else {
    bool foundDiff = false;
    if (src.m_dists[m_dim] != dest->m_dists[m_dim]) {
      DimVec destDims =  dest->m_dists[m_dim].DistEntryDims();
      for (Dim dim = 0; dim < dest->m_numDims; ++dim) {
	if (dim != m_dim) {
	  DistEntry srcDistEntry = src.m_dists[dim];
	  if (dest->m_dists[dim] != srcDistEntry)
	    foundDiff = true;
	  DimVec srcDims = srcDistEntry.DistEntryDims();
	  DimVecIter iter = srcDims.begin();
	  for(; iter != srcDims.end(); ++iter) {
	    DimVecIter iter2 = destDims.begin();
	    for (; iter2 != destDims.end(); ++iter2) {
	      if (*iter == *iter2)
		return false;
	    }
	  }
	}
      }
      return foundDiff;
    }
    else {
      return false;
    }
  }
}

void SplitRedistribs::Apply(Node *node) const
{
  RedistNode *orig = (RedistNode*)node;
  DistType one = orig->InputDataType(0).GetDist();
  const DistType *two = &(orig->m_info.GetDist());

  one.m_dists[m_dim] = two->m_dists[m_dim];


  RedistNode *newRedist = new RedistNode(one, orig->m_info.GetPerm(), 
					 orig->m_align, orig->m_alignModes, orig->m_alignModesSrc);
  newRedist->AddInput(orig->Input(0), orig->InputConnNum(0));
  node->m_poss->AddNode(newRedist);

  if (one == newRedist->InputDataType(0).GetDist())
    throw;


  RedistNode *newRedist2 = new RedistNode(*two, orig->m_info.GetPerm(), 
					 orig->m_align, orig->m_alignModes, orig->m_alignModesSrc);
  newRedist2->AddInput(newRedist, 0);
  node->m_poss->AddNode(newRedist2);

  if (*two == newRedist2->InputDataType(0).GetDist()) {
    cout << newRedist2->GetInputNameStr(0) << endl;
    throw;
  }


  node->RedirectChildren(newRedist2, 0);
  node->m_poss->DeleteChildAndCleanUp(node);

}

bool SingleIndexAllToAll::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
    throw;
  else if (srcType.m_numDims <= m_dim)
    return false;

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = destType.m_dists[m_dim];
  if (srcEntry == destEntry)
    return false;


  if (srcEntry.IsStar() || destEntry.IsStar())
    return false;
  
  DimVec srcDims = srcEntry.DistEntryDims();
  DimVec destDims = destEntry.DistEntryDims();

  if (IsPrefix(srcDims,destDims) || IsPrefix(destDims,srcDims))
    return false;

    DimSet srcSet;
    srcSet.insert(srcDims.begin(), srcDims.end());
    DimSet destSet;
    destSet.insert(destDims.begin(), destDims.end());

  if (srcDims.size() == destDims.size()) {
    //if the src and dest have the same grid modes, then this
    // is already a single-mode AllToAll
    if (includes(srcSet.begin(), srcSet.end(),
		 destSet.begin(), destSet.end()))
      return false;
  }

  //Now check that no other src or dest entry
  //  uses the same grid mode as the src or dest m_dim entry
  DimSet fullSet = srcSet;
  fullSet.insert(destSet.begin(), destSet.end());
  for (Dim dim = 0; dim < srcType.m_numDims; ++dim) {
    if (dim != m_dim) {
      DimSet tmp = srcType.m_dists[dim].DistEntryDimSet();
      DimVec intersection;
      intersection.resize(tmp.size());
      DimVecIter iter = set_intersection(tmp.begin(), tmp.end(),
					 fullSet.begin(), fullSet.end(),
					 intersection.begin());
      if (iter - intersection.begin())
	return false;

      intersection.clear();
      intersection.resize(max(tmp.size(),fullSet.size()));
      tmp = destType.m_dists[dim].DistEntryDimSet();
      iter = set_intersection(tmp.begin(), tmp.end(),
			      fullSet.begin(), fullSet.end(),
			      intersection.begin());
      if (iter - intersection.begin())
	return false;
    }
  }

  return true;
}

void SingleIndexAllToAll::Apply(Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();

  DistType type1 = srcType;
  DistEntry entry1 = type1.m_dists[m_dim];
  DimVec vec1 = entry1.DistEntryDims();
  DimSet set1;
  set1.insert(vec1.begin(), vec1.end());
  
  DimVec destVec = redist->m_info.GetDist().m_dists[m_dim].DistEntryDims();
  DimVecIter iter = destVec.begin();
  for(; iter != destVec.end(); ++iter) {
    Dim dim = *iter;
    if (set1.find(dim) == set1.end()) {
      vec1.push_back(dim);
    }
  }
  type1.m_dists[m_dim].DimsToDistEntry(vec1);

  bool skipFirstRedist = false;
  RedistNode *redist1 = NULL;
  if (type1 == srcType) {
    skipFirstRedist = true;
  }
  else {
    redist1 = new RedistNode(type1, redist->m_info.GetPerm(), 
					 redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
    redist1->AddInput(redist->Input(0), redist->InputConnNum(0));
    node->m_poss->AddNode(redist1);
  }

  DistType type2 = redist->m_info.GetDist();
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
  
  RedistNode *redist2 = new RedistNode(type2, redist->m_info.GetPerm(),
				       redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
  if (skipFirstRedist)
    redist2->AddInput(redist->Input(0), redist->InputConnNum(0));
  else
    redist2->AddInput(redist1, 0);
  node->m_poss->AddNode(redist2);

  if (type2 == redist2->InputDataType(0).GetDist())
    throw;

  RedistNode *redist3 = new RedistNode(redist->m_info.GetDist(), redist->m_info.GetPerm(),
				       redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
  redist3->AddInput(redist2, 0);
  node->m_poss->AddNode(redist3);


  node->RedirectChildren(redist3, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

void PrintVec(DimVec vec)
{
  DimVecIter iter = vec.begin();
  for(; iter != vec.end(); ++iter) {
    cout << *iter << " ";
  }
  cout << endl;
}


bool DoubleIndexAllToAll::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  if (srcType.m_numDims != destType.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_dim)
    return false;

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = destType.m_dists[m_dim];
  if (srcEntry == destEntry)
    return false;
  DimVec srcVec = srcEntry.DistEntryDims();
  DimSet destSet = destEntry.DistEntryDimSet();

  for(Dim dim = m_dim+1; dim < srcType.m_numDims; ++dim) {
    DistEntry srcEntry2 = srcType.m_dists[dim];
    DistEntry destEntry2 = destType.m_dists[dim];
    if (srcEntry2 != destEntry2) {
      DimVec srcVec2 = srcEntry2.DistEntryDims();
      DimSet destSet2 = destEntry2.DistEntryDimSet();
      if (srcVec.empty() && srcVec2.empty())
	continue;
      
      //if at least one grid mode at the end of srcVec is in destSet2
      // and at least one grid mode at the end of srcVec2 is in destSet
      if ((!srcVec.empty() &&  destSet2.find(srcVec[srcVec.size() - 1]) != destSet2.end()) ||
	  (!srcVec2.empty() && destSet.find(srcVec2[srcVec2.size() - 1]) != destSet.end())) {
	{
	  DimVec suff1;
	  DimVecRevIter riter = srcVec.rbegin();
	  for(; riter != srcVec.rend(); ++riter) {
	    Dim d = *riter;
	    if (destSet2.find(d) != destSet2.end()) {
	      suff1.insert(suff1.begin(), d);
	    }
	    else
	      break;
	  }

	  DimVec suff2;
	  riter = srcVec2.rbegin();
	  for(; riter != srcVec2.rend(); ++riter) {
	    Dim d = *riter;
	    if (destSet.find(d) != destSet.end()) {
	      suff2.insert(suff2.begin(), d);
	    }
	    else
	      break;
	  }

	  DistType intType = srcType;

	  int size = suff1.size();
	  if (size) {
	    DimVec int1;
	    if (size < srcVec.size())
	      int1.insert(int1.begin(), srcVec.begin(), srcVec.begin() + srcVec.size() - suff1.size());
	    int1.insert(int1.end(), suff2.begin(), suff2.end());
	    intType.m_dists[m_dim].DimsToDistEntry(int1);
	  }

	  size = suff2.size();
	  if (size) {
	    DimVec int2;
	    if (size < srcVec2.size())
	      int2.insert(int2.begin(), srcVec2.begin(), srcVec2.begin() + srcVec2.size() - suff2.size());
	    int2.insert(int2.end(), suff1.begin(), suff1.end());
	    intType.m_dists[dim].DimsToDistEntry(int2);
	  }

	  if (srcType != intType) {
	    return true;
	  }
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
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  if (srcType.m_numDims != destType.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_dim)
    throw;

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = destType.m_dists[m_dim];
  if (srcEntry == destEntry)
    throw;
  DimVec srcVec = srcEntry.DistEntryDims();
  DimSet destSet = destEntry.DistEntryDimSet();

  for(Dim dim = m_dim+1; dim < srcType.m_numDims; ++dim) {
    DistEntry srcEntry2 = srcType.m_dists[dim];
    DistEntry destEntry2 = destType.m_dists[dim];
    if (srcEntry2 != destEntry2) {
      DimVec srcVec2 = srcEntry2.DistEntryDims();
      DimSet destSet2 = destEntry2.DistEntryDimSet();
      if (srcVec.empty() && srcVec2.empty())
	continue;
      
      //if at least one grid mode at the end of srcVec is in destSet2
      // and at least one grid mode at the end of srcVec2 is in destSet
      if ((!srcVec.empty() &&  destSet2.find(srcVec[srcVec.size() - 1]) != destSet2.end()) ||
	  (!srcVec2.empty() && destSet.find(srcVec2[srcVec2.size() - 1]) != destSet.end())) 
	{
	  DimVec suff1;
	  DimVecRevIter riter = srcVec.rbegin();
	  for(; riter != srcVec.rend(); ++riter) {
	    Dim d = *riter;
	    if (destSet2.find(d) != destSet2.end()) {
	      suff1.insert(suff1.begin(), d);
	    }
	    else
	      break;
	  }

	  DimVec suff2;
	  riter = srcVec2.rbegin();
	  for(; riter != srcVec2.rend(); ++riter) {
	    Dim d = *riter;
	    if (destSet.find(d) != destSet.end()) {
	      suff2.insert(suff2.begin(), d);
	    }
	    else
	      break;
	  }

	  DistType intType = srcType;

	  int size = suff1.size();
	  if (size) {
	    DimVec int1;
	    if (size < srcVec.size())
	      int1.insert(int1.begin(), srcVec.begin(), srcVec.begin() + srcVec.size() - suff1.size());
	    int1.insert(int1.end(), suff2.begin(), suff2.end());
	    intType.m_dists[m_dim].DimsToDistEntry(int1);
	  }

	  size = suff2.size();
	  if (size) {
	    DimVec int2;
	    if (size < srcVec2.size())
	      int2.insert(int2.begin(), srcVec2.begin(), srcVec2.begin() + srcVec2.size() - suff2.size());
	    int2.insert(int2.end(), suff1.begin(), suff1.end());
	    intType.m_dists[dim].DimsToDistEntry(int2);
	  }

	  if (srcType != intType) {
	    RedistNode *intRedist = new RedistNode(intType,
						   redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
	    intRedist->AddInput(redist->Input(0), redist->InputConnNum(0));
		  
	    RedistNode *finalRedist = new RedistNode(destType,
						     redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
	    finalRedist->AddInput(intRedist, 0);

	    redist->RedirectAllChildren(finalRedist);

	    redist->m_poss->AddNode(intRedist);
	    redist->m_poss->AddNode(finalRedist);

	    redist->m_poss->DeleteChildAndCleanUp(redist);


	    return;
	  }
	}
    }
  }

  throw;

}

bool DoubleIndexAllToAllPrefix::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  if (srcType.m_numDims != destType.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_dim)
    return false;

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = destType.m_dists[m_dim];
  if (srcEntry == destEntry)
    return false;
  const DimVec srcVec = srcEntry.DistEntryDims();
  const DimVec destVec = destEntry.DistEntryDims();


  //[(x|_|y),(z|_|w)] <- [(z),(x)]

  for(Dim dim = m_dim+1; dim < srcType.m_numDims; ++dim) {
    DistEntry srcEntry2 = srcType.m_dists[dim];
    DistEntry destEntry2 = destType.m_dists[dim];
    if (srcEntry2 != destEntry2) {
      DimVec srcVec2 = srcEntry2.DistEntryDims();
      DimVec destVec2 = destEntry2.DistEntryDims();
      if (IsPrefix(srcVec, destVec2) && IsPrefix(srcVec2, destVec))
	return true;
      }
    }
  return false;
}

void DoubleIndexAllToAllPrefix::Apply(Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  //  cout << srcType.PrettyStr() << " -> " << redist->DataType(0).m_dist.PrettyStr() << endl;

  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
    throw;
  else if (srcType.m_numDims <= m_dim)
    throw;

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = redist->m_info.GetDist().m_dists[m_dim];
  if (srcEntry == destEntry)
    throw;
  const DimVec srcVec = srcEntry.DistEntryDims();
  const DimVec destVec = destEntry.DistEntryDims();

  for(Dim dim = m_dim+1; dim < srcType.m_numDims; ++dim) {
    DistEntry srcEntry2 = srcType.m_dists[dim];
    DistEntry destEntry2 = redist->m_info.GetDist().m_dists[dim];
    if (srcEntry2 != destEntry2) {
      DimVec srcVec2 = srcEntry2.DistEntryDims();
      DimVec destVec2 = destEntry2.DistEntryDims();
      if (IsPrefix(srcVec, destVec2) && IsPrefix(srcVec2, destVec)) {
	//[(x|_|y),(z|_|w)] <- [(z),(x)]
	// to
	//[(x),(z)] <- [(z),(x)]
	//[(x|_|y),(z|_|w)] <- [(x),(z)]

	DistType intType = srcType;
	intType.m_dists[m_dim] = srcEntry2;
	intType.m_dists[dim] = srcEntry;
	
	RedistNode *newRedist = new RedistNode(intType, redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
	newRedist->AddInput(redist->Input(0), redist->InputConnNum(0));

	redist->m_poss->AddNode(newRedist);

	redist->ChangeInput2Way(redist->Input(0), redist->InputConnNum(0), newRedist, 0);

	return;
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
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
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
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = destType.m_dists[m_dim];

  DimVec destDims = destEntry.DistEntryDims();
  DimVec srcDims = srcEntry.DistEntryDims();

  srcDims.pop_back();

  DistType intType = srcType;
  intType.m_dists[m_dim].DimsToDistEntry(srcDims);

  RedistNode *redist1 = new RedistNode(intType, redist->m_info.GetPerm(),
				       redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
  redist1->AddInput(redist->Input(0), redist->InputConnNum(0));
  node->m_poss->AddNode(redist1);

  if (intType == redist1->InputDataType(0).GetDist())
    throw;

  RedistNode *redist2 = new RedistNode(destType, redist->m_info.GetPerm(),
				       redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
  redist2->AddInput(redist1, 0);
  node->m_poss->AddNode(redist2);

  if (destType == redist2->InputDataType(0).GetDist())
    throw;

  node->RedirectChildren(redist2, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

bool SplitAllAllGathers::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
    throw;

  int allGathers = 0;
  bool foundOther = false;

  for (Dim dim = 0; dim < srcType.m_numDims; ++dim) {
    const DistEntry &src = srcType.m_dists[dim];
    const DistEntry &dest = destType.m_dists[dim];
    if (src != dest) {
      if (dest.IsStar()) {
	++allGathers;
      }
      else if (IsPrefix(dest.DistEntryDims(),
			src.DistEntryDims())) 
	{
	  ++allGathers;
	}
      else
	foundOther = true;
    }
  }
  return (allGathers > 1) && foundOther;
}

void SplitAllAllGathers::Apply(Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  DistType intType = srcType;

  for (Dim dim = 0; dim < srcType.m_numDims; ++dim) {
    const DistEntry &src = srcType.m_dists[dim];
    const DistEntry &dest = destType.m_dists[dim];
    if (src != dest) {
      if (dest.IsStar()) {
	intType.m_dists[dim] = dest;
      }
      else if (IsPrefix(dest.DistEntryDims(),
			src.DistEntryDims())) 
	{
	  intType.m_dists[dim] = dest;
	}
    }
  }

  RedistNode *redist1 = new RedistNode(intType, redist->m_info.GetPerm(),
				       redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
  redist1->AddInput(redist->Input(0), redist->InputConnNum(0));
  node->m_poss->AddNode(redist1);

  if (intType == redist1->InputDataType(0).GetDist())
    throw;

  RedistNode *redist2 = new RedistNode(destType, redist->m_info.GetPerm(),
				       redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
  redist2->AddInput(redist1, 0);
  node->m_poss->AddNode(redist2);

  if (destType == intType)
    throw;

  node->RedirectChildren(redist2, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

bool CombineDisappearingModes::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
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
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
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

  RedistNode *newRedist = new RedistNode(intType, redist->m_info.GetPerm(),
                                         redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
  newRedist->AddInput(redist->Input(0), redist->InputConnNum(0));
  redist->m_poss->AddNode(newRedist);

  if (intType == newRedist->InputDataType(0).GetDist())
    throw;
  
  redist->ChangeInput2Way(redist->Input(0), redist->InputConnNum(0),
			  newRedist, 0);    
}

bool PermuteDistribution::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;

  const RedistNode *redist = (RedistNode*)node;

  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
    throw;
  else if (srcType.m_numDims <= m_srcDim || srcType.m_numDims <= m_destDim)
    return false;


  const DistEntry &srcEntry = srcType.m_dists[m_srcDim];
  const DistEntry &destEntry = destType.m_dists[m_destDim];

  /*
  if (m_srcDim == 1 && m_destDim == 0) {
    if (srcEntry.m_val == 8 && destEntry.m_val == 6 && srcType.m_dists[0].m_val == 11) {
      cout << srcEntry.PrettyStr() << endl;
      cout << destEntry.PrettyStr() << endl;
      cout << srcType.PrettyStr() << endl;
      cout << destType.PrettyStr() << endl;
    }
  }
  */


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

  /*  
  cout << srcType.PrettyStr() << endl;
  cout << type.PrettyStr() << endl;
  cout << destType.PrettyStr() << endl;
  */

  return DistTypeNotEqual(type, srcType);
}

void PermuteDistribution::Apply(Node *node) const
{  
  // Didn't know if this was actually needed
  //check / propagate m_permutation as needed
  //  throw;
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  RedistNode *redist = (RedistNode*)node;

  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
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

  RedistNode *newRedist = new RedistNode(type,
                                         redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
  newRedist->AddInput(redist->Input(0), redist->InputConnNum(0));
  redist->m_poss->AddNode(newRedist);
  
  RedistNode *newRedist2 = new RedistNode(redist->m_info.GetDist(),
                                          redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
  newRedist2->AddInput(newRedist,0);
  redist->m_poss->AddNode(newRedist2);
  
  redist->RedirectAllChildren(newRedist2);
  
  redist->m_poss->DeleteChildAndCleanUp(redist);
}

bool GetAllToAllPattern(const DistType &srcType, 
			const DistType &destType,
			DistEntryVec *gridModesInvolved)
{
  map<Dim,DistEntry> srcModes;
  map<DistEntry,Dim,DistEntryCompare> destModes;
  typedef pair<Dim,DistEntry> SrcModePair;
  typedef pair<DistEntry,Dim> DestModePair;

  Dim numDims = srcType.m_numDims;
  
  EntrySet srcSuffixes, destSuffixes;
  for (Dim dim = 0; dim < numDims; ++dim) {
    DistEntry srcDistEntry = srcType.m_dists[dim];
    DistEntry destDistEntry = destType.m_dists[dim];
    DimVec srcDims = srcDistEntry.DistEntryDims();
    DimVec destDims = destDistEntry.DistEntryDims();
    DimVec suff1, suff2;
    if (srcDistEntry != destDistEntry) {
      GetDifferentSuffix(srcDims, destDims, suff1, suff2);
      if (suff1.empty() || suff2.empty()) {
	return false;
      }
      DistEntry srcSuff;
      srcSuff.DimsToDistEntry(suff1);
      srcSuffixes.insert(srcSuff);
      DistEntry destSuff;
      destSuff.DimsToDistEntry(suff2);
      destSuffixes.insert(destSuff);
    }
  }

  EntrySetIter srcIter = srcSuffixes.begin();
  for(; srcIter != srcSuffixes.end(); ++srcIter) {
    DistEntry entry = *srcIter;
    if (!destSuffixes.erase(entry))
      return false;
  }
  if (!destSuffixes.empty())
    return false;

  if (gridModesInvolved) {
    map<Dim,DistEntry>::iterator srcModeIter = srcModes.begin();
    for(; srcModeIter != srcModes.end(); ++srcModeIter) {
      gridModesInvolved->push_back(srcModeIter->second);
    }
  }

  return true;
}

bool GetAllGatherPattern(const DistType &srcType,
			 const DistType &destType,
			 DistEntryVec *gridModes)
{
  const DimSet destUsedDims = destType.UsedGridDims();

  Dim numDims = srcType.m_numDims;

  if (numDims != destType.m_numDims)
    return false;

  int count = 0;
  
  for (Dim dim = 0; dim < numDims; ++dim) {
    DistEntry srcDistEntry = srcType.m_dists[dim];
    DistEntry destDistEntry = destType.m_dists[dim];
    if (srcDistEntry != destDistEntry) {
      DimVec srcDims = srcDistEntry.DistEntryDims();
      DimVec destDims = destDistEntry.DistEntryDims();
      DimVec suff;
      if (srcDims.size() > destDims.size()) {
	GetSuffix(destDims, srcDims, suff);
	if (suff.size() + destDims.size() != srcDims.size())
	  return false;
	DimVecIter iter = suff.begin();
	for(; iter != suff.end(); ++iter) {
	  if (destUsedDims.find(*iter) != destUsedDims.end())
	    return false;
	}
	if (gridModes) {
	  DistEntry entry;
	  entry.DimsToDistEntry(suff);
	  gridModes->push_back(entry);
	}
	++count;
      }
      else
	return false;
    }
  }

#if ALLMULTIMODEALLGATHER
  return count;
#else
  return count==1;
#endif
}

bool GetLocalRedistPattern(const DistType &srcType,
			   const DistType &destType,
			   DimVec *indices,
			   DistEntryVec *gridModes)
{
  const DimSet srcUsedDims = srcType.UsedGridDims();

  Dim numDims = srcType.m_numDims;

  if (numDims != destType.m_numDims)
    return false;

  int count = 0;
  
  for (Dim dim = 0; dim < numDims; ++dim) {
    DistEntry srcDistEntry = srcType.m_dists[dim];
    DistEntry destDistEntry = destType.m_dists[dim];
    if (srcDistEntry != destDistEntry) {
      DimVec srcDims = srcDistEntry.DistEntryDims();
      DimVec destDims = destDistEntry.DistEntryDims();
      DimVec suff;
      if (destDims.size() > srcDims.size()) {
	GetSuffix(srcDims, destDims, suff);
	if (suff.size() + srcDims.size() != destDims.size())
	  return false;
	DimVecIter iter = suff.begin();
	for(; iter != suff.end(); ++iter) {
	  if (srcUsedDims.find(*iter) != srcUsedDims.end())
	  return false;
	}
	if (indices) {
	  indices->push_back(dim);
	  DistEntry entry;
	  entry.DimsToDistEntry(suff);
	  gridModes->push_back(entry);
	}
	++count;
      }
      else
	return false;
    }
  }

  return count != 0;
}

#endif



