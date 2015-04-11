/*
 This file is part of DxTer.
 DxTer is a prototype using the Design by Transformation (DxT)
 approach to program generation.
 
 Copyright (C) 2015, The University of Texas and Bryan Marker
 
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
#include "tempVarNode.h"

#if DOTENSORS

#include "tensorSumScatter.h"

extern bool M_allowSquareGridOpt;


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
                        DimVec *gridModesInvolved);

bool GetAllGatherPattern(const DistType &srcType,
                         const DistType &destType,
                         DimVec *gridModes);

bool GetLocalRedistPattern(const DistType &srcType,
                           const DistType &destType);

bool GetPermPattern(const DistType &srcType,
		    const DistType &destType,
		    DimVec *gridModes);


RedistNode::RedistNode()
{
  m_info.SetToDefault(0);
}

RedistNode::RedistNode(const DistType &destType, const string &align,
                       const DimVec &alignModes, const DimVec &alignModesSrc)
{
  m_info.SetDistAndClearPerm(destType);
  m_align = align;
  m_alignModes = alignModes;
  m_alignModesSrc = alignModesSrc;
}


RedistNode::RedistNode(const DistType &destType)
{
  m_info.SetDistAndClearPerm(destType);
  m_align = "NONE";
}

RedistNode::RedistNode(const DistType &destType, const Permutation &perm, const string &align,
                       const DimVec &alignModes, const DimVec &alignModesSrc)
: RedistNode(destType, align, alignModes, alignModesSrc)
{
  m_info.SetPerm(perm);
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
    return (string)"RedistNode to " +  m_info.GetDist().QuickStr() 
      + m_info.GetPerm().Str() + m_align;
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
      m_name = (string)"RedistNode to " +  m_info.GetDist().QuickStr() + m_info.GetPerm().Str() + m_align;
    DLANode *parent = (DLANode*)Input(0);
    parent->Prop();
    
    
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
    
    DimVec gridModesInvolved;
    
    if (GetAllGatherPattern(m_srcType, m_info.GetDist(),
                            &gridModesInvolved)) {
      unsigned int numProcs = 1;
      DimVecIter gridModeIter = gridModesInvolved.begin();
      for(; gridModeIter != gridModesInvolved.end(); ++gridModeIter) {
        numProcs *= GridLens[*gridModeIter];
      }

      const unsigned int totNumIters = m_lsizes[0]->NumSizes();
      for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
        Cost temp = 1;
        for (Dim dim = 0; dim < numDims; ++dim) {
          temp *= (*(m_lsizes[dim]))[iteration];
        }

        m_cost += AllGather(temp, numProcs);
        m_cost += (PSIR+PSIW)*(temp + temp / numProcs);
      }

      return;
    }
    
    DimVec indices;
    if (GetLocalRedistPattern(m_srcType, m_info.GetDist())) {
      //Local mem copy
      m_cost = (PSIR+PSIW)*(TotalNumberOfLocalElements(0));
      return;
    }

    if (GetPermPattern(m_srcType, m_info.GetDist(),
		       &gridModesInvolved)) {
      Cost comm = 0;
      Cost mov = 0;
      const unsigned int totNumIters = m_lsizes[0]->NumSizes();
      for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
        Cost tempOut = 1;
        for (Dim dim = 0; dim < numDims; ++dim) {
          tempOut *= (*(m_lsizes[dim]))[iteration];
        }
        comm += SendRecv(tempOut);
	mov += (PSIR+PSIW)*(2*tempOut);
      }
      m_cost = comm + mov;
      //m_cost = comm;
      return;
    }    
    
    if (GetAllToAllPattern(m_srcType, m_info.GetDist(),
                           &gridModesInvolved)) {
      unsigned int numProcs = 1;
      DimVecIter gridModeIter = gridModesInvolved.begin();
      for(; gridModeIter != gridModesInvolved.end(); ++gridModeIter) {
        numProcs *= GridLens[*gridModeIter];
      }
      const unsigned int totNumIters = m_lsizes[0]->NumSizes();
      Cost comm = 0;
      Cost mov = 0;
      for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
        Cost tempOut = 1;
        for (Dim dim = 0; dim < numDims; ++dim) {
          tempOut *= (*(m_lsizes[dim]))[iteration];
        }
        comm += AllToAll(tempOut, numProcs);
	mov += (PSIR+PSIW)*(2*tempOut);
      }
      m_cost = comm + mov;
      //m_cost = comm;

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
      
      const unsigned int totNumIters = m_lsizes[0]->NumSizes();
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
        Cost tempOut = 1;
	Cost tempIn = 1;
        for (Dim dim = 0; dim < numDims; ++dim) {
          tempOut *= (*(m_lsizes[dim]))[iteration];
	  tempIn *= (*InputLocalLen(0, dim))[iteration];
        }
        m_cost += AllToAll(tempOut, numProcs);
        m_cost += (PSIR+PSIW)*(tempIn + tempOut);
      }
    }
    else
      m_cost = 0;
    
  }
}

bool RedistNode::HasSamePerm(const RedistNode *redist) const
{
  if (m_info.GetPerm() != redist->m_info.GetPerm())
    {
      return false;
    }
  else
    return true;

}

bool RedistNode::HasSamePermAndAlign(const RedistNode *redist) const
{
  if (!HasSamePerm(redist)
      || !HasSameAlign(redist))
    {
      return false;
      
    }
  else
    return true;

}


bool RedistNode::HasSameAlign(const RedistNode *redist) const
{
  if (m_align != redist->m_align
      || m_alignModes != redist->m_alignModes
      || m_alignModesSrc != redist->m_alignModesSrc)
    {
      return false;
      
    }
  else
    return true;

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

  if (GetPermPattern(m_srcType, m_info.GetDist(),
			  NULL))
    return true;
  
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
  
  if (GetLocalRedistPattern(m_srcType, m_info.GetDist()))
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

const SizeList* RedistNode::Len(ConnNum num, Dim dim) const
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

const SizeList* RedistNode::LocalLen(ConnNum num, Dim dim) const
{
  if (num > 0)
    throw;
  //Input perm is accounted for in building the cache
  //  if (InputDataType(0).HasPerm())
  //    throw;
  if (!m_isArray) {
    return m_lsizes[0];
  }
  else {
    if (m_info.HasPerm())
      return m_lsizes[m_info.GetPerm().MapFinishToStart(dim)];
    else
      return m_lsizes[dim];
  }
}

void RedistNode::ClearDataTypeCache()
{
  m_lsizes.clear();
}

void RedistNode::BuildDataTypeCache()
{
  if (!m_lsizes.empty())
    return;
  
  //If the input has a permutation, 
  //it has to be undone to get local sizes of inputs

  DLANode *in = (DLANode*)Input(0);
  ConnNum num = InputConnNum(0);
  Dim numDims = in->NumDims(num);

  Permutation perm = in->DataType(num).GetPerm();

  if (numDims) {
    m_isArray = true;
    for (Dim dim = 0; dim < numDims; ++dim) {
      m_lsizes.push_back(GetLocalSizes(in->Len(num,perm.MapStartToFinish(dim)), m_info.GetDist().m_dists[dim]));
    }
  }
  else {
    m_isArray = false;
    m_lsizes.push_back(in->Len(num,0));
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
  
  if (!HasNoAlign()) {
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
  
  DimVec gridModesInvolved;
  
  
  DimVec indices;
  if (GetAllGatherPattern(m_srcType, m_info.GetDist(),
                          &gridModesInvolved)) {
    if (gridModesInvolved.empty())
      throw;
    *out << outName << ".AllGatherRedistFrom( "
    << inName << ", "
    << ModeArrayVarName(gridModesInvolved) << " );\n";

#if 0
      double numProcs = 1;
      DimVecIter gridModeIter = gridModesInvolved.begin();
      for(; gridModeIter != gridModesInvolved.end(); ++gridModeIter) {
        numProcs *= GridLens[*gridModeIter];
      }
      *out << "for previous one\n";
      *out << "numProcs " << numProcs << endl;

      Size size = 1;
      const unsigned int totNumIters = m_lsizes[0]->NumSizes();
      for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
        Cost temp = 1;
        for (Dim dim = 0; dim < numDims; ++dim) {
          size *= (*(m_lsizes[dim]))[iteration];
        }
      }
      *out << "on " << size << endl;
      *out << log2(numProcs) * ALPHA << " alpha\n";
      *out << ((numProcs - ONE) / numProcs) << endl;
      *out << ((numProcs - ONE) / numProcs) * size * BETA;
      *out << " beta\n";
#endif // 0

    return;
  }
  
  if (GetLocalRedistPattern(m_srcType, m_info.GetDist())) {
    *out << outName << ".LocalRedistFrom( "
    << inName << " );\n";
    return;
  }

  if (GetPermPattern(m_srcType, m_info.GetDist(),
                         &gridModesInvolved)) {
    *out << outName << ".PermutationRedistFrom( "
	 << inName << ", "
	 << ModeArrayVarName(gridModesInvolved) << " );\n";
    return;
  }
  
  if (GetAllToAllPattern(m_srcType, m_info.GetDist(),
                         &gridModesInvolved)) {
    if (gridModesInvolved.empty()) {
      GetAllToAllPattern(m_srcType, m_info.GetDist(),
                         &gridModesInvolved);
      throw;
    }
    *out << outName << ".AllToAllRedistFrom( "
    << inName << ", "
    << ModeArrayVarName(gridModesInvolved) << " );\n";
#if 0
    *out << "for that previous one:\n";
    double numProcs = 1;
    DimVecIter gridModeIter = gridModesInvolved.begin();
    for(; gridModeIter != gridModesInvolved.end(); ++gridModeIter) {
      numProcs *= GridLens[*gridModeIter];
    }
    *out << "numProcs = " << numProcs << endl;
    Size size = 0;
    double cost = 0;
      const unsigned int totNumIters = m_lsizes[0]->NumSizes();
      for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
        Cost tempOut = 1;
	Cost tempIn = 1;
        for (Dim dim = 0; dim < numDims; ++dim) {
          tempOut *= (*(m_lsizes[dim]))[iteration];
	  tempIn *= (*InputLocalLen(0, dim))[iteration];
        }
	size += tempOut;
	//        cost += AllToAll(tempOut, numProcs);
	//       m_cost += (PSIR+PSIW)*(tempIn + tempOut);
      }
      *out << "size " << size << endl;
      *out << (numProcs-1) * ALPHA << " alpha\n";
      *out << ((numProcs-1)/numProcs) << endl;
      *out << ((numProcs-1)/numProcs) * (size * BETA);
      *out << " beta\n";
#endif // 0
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
    << inName << ", "
    << ModeArrayVarName(src) << " );\n";
    return;
  }
  else
    throw;
}

bool RedistNode::HasNoAlign() const
{
  if (m_alignModes.empty() || m_align == "NONE")
    return true; // to superscede any of the following
  else if (m_align != GetInputNameStr(0))
    return false;
  else {
    //aligning to input, so make sure the alignment isn't the identity
    if (m_alignModes != m_alignModesSrc)
      return false;
    else
      return true;
  }
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
  
  DimVec gridModesInvolved;
  
  if (GetAllGatherPattern(m_srcType, m_info.GetDist(),
                          &gridModesInvolved)) {
    {
      Var var(ModeArrayVarType, gridModesInvolved);
      set.insert(var);
    }
    return;
  }
  
  DimVec indices;
  if (GetLocalRedistPattern(m_srcType, m_info.GetDist())) {
    return;
  }
  
  if (GetPermPattern(m_srcType, m_info.GetDist(),
                         &gridModesInvolved)) {
    {
      Var var(ModeArrayVarType, gridModesInvolved);
      set.insert(var);
    }
    
    return;
  }


  if (GetAllToAllPattern(m_srcType, m_info.GetDist(),
                         &gridModesInvolved)) {
    {
      Var var(ModeArrayVarType, gridModesInvolved);
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
  bool onlyOneChild = true;
  while(redistNode->Input(0)
        && (redistNode->Input(0)->GetNodeClass() == RedistNode::GetClass()))
  {
    if (redistNode->m_children.size() > 1)
      onlyOneChild = false;
    redistNode = (RedistNode*)redistNode->Input(0);
    if (DistTypeEqual(redistNode->m_info.GetDist(),*type) 
	&& redistNode->HasSamePerm((RedistNode*)node)) 
      {
	if (!onlyOneChild) {
	  return redistNode->HasSameAlign((RedistNode*)node);
	}
	else {
	  return redistNode->HasNoAlign() || redistNode->HasSameAlign((RedistNode*)node);
	}
      }
    for(ConnNum i = 0; i < redistNode->m_children.size(); ++i) {
      Node *tmp = redistNode->Child(i);
      if (tmp != node && tmp->GetNodeClass() == RedistNode::GetClass()) {
        RedistNode *redist = ((RedistNode*)tmp);
        if (DistTypeEqual(redist->m_info.GetDist(), *type)) {
	  return redistNode->HasSamePermAndAlign((RedistNode*)node);
	}
      }
    }
  }
  if (redistNode->Input(0)
      && DistTypeEqual(((DLANode*)(redistNode->Input(0)))->DataType(redistNode->InputConnNum(0)).GetDist(),*type)
      && ((RedistNode*)node)->HasNoAlign()
      && !((RedistNode*)node)->m_info.HasPerm())
    return true;
  return false;
}

void RemoveWastedRedist::Apply(Node *node) const
{
  RedistNode *redistNode = (RedistNode*)node;
  const DistType* type = &(redistNode->m_info.GetDist());
  bool onlyOneChild = true;
  while(redistNode->Input(0)
        && (redistNode->Input(0)->GetNodeClass() == RedistNode::GetClass()))
  {
    if (redistNode->m_children.size() > 1)
      onlyOneChild = false;
    redistNode = (RedistNode*)redistNode->Input(0);
    if (DistTypeEqual(redistNode->m_info.GetDist(), *type)
	&& redistNode->HasSamePerm((RedistNode*)node))
    {
      node->RedirectChildren(redistNode, 0);
      if (!redistNode->HasSameAlign((RedistNode*)node)) {
	if (!onlyOneChild)
	  throw;
	else {
	  redistNode->m_align = ((RedistNode*)node)->m_align;
	  redistNode->m_alignModes = ((RedistNode*)node)->m_alignModes;
	  redistNode->m_alignModesSrc = ((RedistNode*)node)->m_alignModesSrc;	  
	}
      }
      node->m_poss->DeleteChildAndCleanUp(node);
      return;
    }
    for(ConnNum i = 0; i < redistNode->m_children.size(); ++i) {
      Node *tmp = redistNode->Child(i);
      if (tmp != node && tmp->GetNodeClass() == RedistNode::GetClass()) {
        RedistNode *redist = ((RedistNode*)tmp);
        if (DistTypeEqual(redist->m_info.GetDist(),*type)) {
	  node->RedirectChildren(tmp, 0);
	  node->m_poss->DeleteChildAndCleanUp(node);
	  return;
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

bool CombinePermuteRedists::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  if (node->m_children.size() != 1)
    return false;
  const DistType &srcType = redist->InputDataType(0).GetDist();

  if (node->Child(0)->GetNodeClass() != RedistNode::GetClass())
    return false;

  const RedistNode *redist2 = (RedistNode*)(node->Child(0));
  const DistType &destType2 = redist2->m_info.GetDist();
  if (!redist->HasNoAlign() && !redist->HasSameAlign(redist2))
    return false;
  if (!GetPermPattern(srcType, destType2, NULL))
    return false;
  else
    return true;
  
}

void CombinePermuteRedists::Apply(Node *node) const
{
  RedistNode *redist = (RedistNode*)node;
  RedistNode *redist2 = (RedistNode*)(node->Child(0));
  redist2->ChangeInput2Way(redist, 0, 
			   redist->Input(0), redist->InputConnNum(0));
  redist->m_poss->DeleteChildAndCleanUp(redist);
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
	&& redist->HasSamePermAndAlign((RedistNode*)output))
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
	&& redist->HasSamePermAndAlign((RedistNode*)output))
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

  if (GetPermPattern(src, *dest, NULL))
    return false;      
  
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
  
  RedistNode *newRedist = new RedistNode(one);
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

  if (GetPermPattern(srcType, destType, NULL))
    return false;
  
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
    redist1 = new RedistNode(type1);
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
  
  RedistNode *redist2 = new RedistNode(type2);
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

  if (GetPermPattern(srcType, destType, NULL))
    return false;
  
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
          
          unsigned int size = suff1.size();
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
        
        unsigned int size = suff1.size();
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
          RedistNode *intRedist = new RedistNode(intType);
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

  if (GetPermPattern(srcType, destType, NULL))
    return false;
  
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

        RedistNode *newRedist = new RedistNode(intType);
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
#if ALLMULTIMODEALLGATHER
    return !GetAllGatherPattern(srcType, destType, NULL);
#else
  return true;
#endif
  
  DimVec destDims = destEntry.DistEntryDims();
  
  if (srcDims.size() < destDims.size())
    return false;
  
  if ((srcDims.size() - destDims.size()) > 1 &&  IsPrefix(destDims, srcDims)) {
#if ALLMULTIMODEALLGATHER
    return !GetAllGatherPattern(srcType, destType, NULL);
#else
    return true;
#endif
  }
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

  RedistNode *redist1 = new RedistNode(intType);
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

  RedistNode *redist1 = new RedistNode(intType);
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

#if ALLMULTIMODEALLGATHER
bool CombineAllGathers::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();
  
  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
    throw;
  
  if (!GetAllGatherPattern(srcType, destType, NULL))
    return false;
  
  if (redist->m_children.size() != 1) {
    return false;
  }
  
  const Node *child = redist->Child(0);
  if (child->GetNodeClass() != RedistNode::GetClass())
    return false;
  
  const RedistNode *childRedist = (RedistNode*)child;

  //  if (!redist->HasSamePerm(childRedist))
  //    return false;

  if (!redist->HasNoAlign() && !redist->HasSameAlign(childRedist))
    return false;
  
  return GetAllGatherPattern(srcType, childRedist->m_info.GetDist(), NULL);
}

void CombineAllGathers::Apply(Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  
  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
    throw;
  
  if (redist->m_children.size() != 1) {
    throw;
  }
  
  Node *child = redist->Child(0);
  if (child->GetNodeClass() != RedistNode::GetClass())
    throw;
  
  RedistNode *childRedist = (RedistNode*)child;
  
  RedistNode *newRedist = new RedistNode(childRedist->m_info.GetDist(),
                                         childRedist->m_info.GetPerm(),
                                         childRedist->m_align,
                                         childRedist->m_alignModes,
                                         childRedist->m_alignModesSrc);
  newRedist->AddInput(redist->Input(0), redist->InputConnNum(0));
  node->RedirectChildren(newRedist, 0);
  node->m_poss->AddNode(newRedist);
  
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif //ALLMULTIMODEALLGATHER

bool CombineDisappearingModes::CanApply(const Node *node) const
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

  if (destType.m_dists[m_dim] == srcType.m_dists[m_dim] ||
      destType.m_dists[m_dim].IsStar() ||
      srcType.m_dists[m_dim].IsStar())
    return false;
  
  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = destType.m_dists[m_dim];

  DimVec srcDims = srcEntry.DistEntryDims();
  DimSet destSet = destEntry.DistEntryDimSet();

  DimVec prefix, suffix;

  DimVecIter iter = srcDims.begin();
  for( ; iter != srcDims.end(); ++iter) {
    if (destSet.find(*iter) == destSet.end()) {
      suffix.push_back(*iter);
    }
    else {
      prefix.push_back(*iter);
    }
  }

  prefix.insert(prefix.end(),
		suffix.begin(), suffix.end());

  return (prefix != srcDims);
    
}

void CombineDisappearingModes::Apply(Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;

  RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();
  
  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = destType.m_dists[m_dim];

  DimVec srcDims = srcEntry.DistEntryDims();
  DimSet destSet = destEntry.DistEntryDimSet();


  DimVec prefix, suffix;

  DimVecIter iter = srcDims.begin();
  for( ; iter != srcDims.end(); ++iter) {
    if (destSet.find(*iter) == destSet.end()) {
      suffix.push_back(*iter);
    }
    else {
      prefix.push_back(*iter);
    }
  }

  prefix.insert(prefix.end(),
		suffix.begin(), suffix.end());
  
  DistType intType = srcType;
  intType.m_dists[m_dim].DimsToDistEntry(prefix);

  RedistNode *newRedist = new RedistNode(intType);
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

  if (GetPermPattern(srcType, destType, NULL))
    return false;
  
  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
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
  
  return DistTypeNotEqual(type, srcType) && DistTypeNotEqual(type,destType);
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

  RedistNode *newRedist = new RedistNode(type);
  newRedist->AddInput(redist->Input(0), redist->InputConnNum(0));
  redist->m_poss->AddNode(newRedist);
  
  RedistNode *newRedist2 = new RedistNode(redist->m_info.GetDist(), redist->m_info.GetPerm(),
                                          redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
  newRedist2->AddInput(newRedist,0);
  redist->m_poss->AddNode(newRedist2);
  
  redist->RedirectAllChildren(newRedist2);
  
  redist->m_poss->DeleteChildAndCleanUp(redist);
}

bool GetAllToAllPattern(const DistType &srcType,
                        const DistType &destType,
                        DimVec *gridModesInvolved)
{
  // For each tensor mode, check the different
  // suffixes in the src and dest.
  // [0123,4567] -> [0276,453]
  // AllToAll with 1 2 3 5 6 7
  // Have basic check to not include [01,23] -> [0,23] (just AllGather)
  // or [0,1] -> [02,1] (local redist)
  // Also, don't include multi-tensor-mode AllToAll when single mode would work
  //   e.g., [0123,456] -> [0321, 654]
  
  const DimSet srcSet = srcType.UsedGridDims();
  const DimSet destSet = destType.UsedGridDims();

  DimSet gridSet;
  
  Dim numDims = srcType.m_numDims;
  
  bool foundSingleModeAllToAll = false;
  bool foundMultiModeAllToAll = false;
  
  for (Dim dim = 0; dim < numDims; ++dim) {
    DistEntry srcEntry = srcType.m_dists[dim];
    DistEntry destEntry = destType.m_dists[dim];
    DimVec srcDims = srcEntry.DistEntryDims();
    DimVec destDims = destEntry.DistEntryDims();
    if (srcEntry != destEntry) {
      DimVec suff1, suff2;
      GetDifferentSuffix(srcDims, destDims, suff1, suff2);
      if (gridModesInvolved) {
        gridSet.insert(suff1.begin(), suff1.end());
        gridSet.insert(suff2.begin(), suff2.end());
      }
      if (suff1.empty() && !suff2.empty()) {
        //Check if this is just a local redist
        // If a grid mode in the dest suffix is found in
        // another src tensor mode, then this isn't
        bool foundAllToAll = false;
        DimVecIter iter = suff2.begin();
        for(; !foundAllToAll && iter != suff2.end(); ++iter) {
          if (srcSet.find(*iter) != srcSet.end()) {
            foundAllToAll = true;
          }
        }
        //just a localredist
        if (!foundAllToAll)
          return false;
        //make sure that the final grid mode is found in another tensor
        // mode. otherwise, a local redist should be split off
        if (srcSet.find(suff2.back()) == srcSet.end()) {
          return false;
        }
        foundMultiModeAllToAll = true;
      }
      else if (!suff1.empty() && suff2.empty()) {
        bool foundAllToAll = false;
        //make sure this isn's just an allgather
        DimVecIter iter = suff1.begin();
        for(; !foundAllToAll && iter != suff1.end(); ++iter) {
          if (destSet.find(*iter) != destSet.end())
            foundAllToAll = true;
        }
        if (!foundAllToAll)
          return false;
        foundMultiModeAllToAll = true;
      }
      else if (!suff1.empty() && !suff2.empty()) {
        //if the modes in suff1 aren't found in different
        // tensor modes, this is a single mode all to all
        bool foundMultiMode = false;
        DimVecIter iter = suff1.begin();
        for(; !foundMultiMode && iter != suff1.end(); ++iter) {
          if (!FoundInDimVec(suff2, *iter)) {
            if (destSet.find(*iter) != destSet.end()) {
              foundMultiMode = true;
            }
          }
        }

        iter = suff2.begin();
        for(; !foundMultiMode && iter != suff2.end(); ++iter) {
          if (!FoundInDimVec(suff1, *iter)) {
            if (srcSet.find(*iter) != srcSet.end()) {
              foundMultiMode = true;
            }
          }
        }

        if (foundMultiMode)
          foundMultiModeAllToAll = true;
        else {
          //this is a single-mode AllToAll
          // if there are multiple single-mode AllToAll,
          // they should be split apart
          if (foundSingleModeAllToAll)
            return false;
          foundSingleModeAllToAll = true;
        }
      }
    }
  }

  if (gridModesInvolved) {
    gridModesInvolved->insert(gridModesInvolved->begin(),
			      gridSet.begin(), gridSet.end());
  }
  
  if (foundSingleModeAllToAll) {
    //If there's a single mode AND multi-mode AllToAll,
    // the single AllToAll should be separated
    if (foundMultiModeAllToAll)
      return false;
    return true;
  }
  else if (foundMultiModeAllToAll)
    return true;
  else
    return false;
}

bool GetAllGatherPattern(const DistType &srcType,
                         const DistType &destType,
                         DimVec *gridModes)
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
          gridModes->insert(gridModes->end(), suff.begin(), suff.end());
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
                           const DistType &destType)
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
        ++count;
      }
      else
        return false;
    }
  }
  
  return count != 0;
}



bool DoubleIndexAllToAll2::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();

  if (GetPermPattern(srcType, destType, NULL))
    return false;
  
  if (srcType.m_numDims != redist->m_info.GetDist().m_numDims)
    throw;
  else if (srcType.m_numDims <= m_dim)
    return false;
  
  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = destType.m_dists[m_dim];
  if (srcEntry == destEntry)
    return false;
  
  for(Dim dim = m_dim+1; dim < srcType.m_numDims; ++dim) {
    DistEntry srcEntry2 = srcType.m_dists[dim];
    DistEntry destEntry2 = destType.m_dists[dim];
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
          if (DistTypeNotEqual(intType, redist->m_info.GetDist()))
            return true;
        }
      }
    }
  }
  return false;
}

void DoubleIndexAllToAll2::Apply(Node *node) const
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
  
  for(Dim dim = m_dim+1; dim < srcType.m_numDims; ++dim) {
    DistEntry srcEntry2 = srcType.m_dists[dim];
    DistEntry destEntry2 = redist->m_info.GetDist().m_dists[dim];
    if (srcEntry2 != destEntry2) {
      DimVec srcVec = srcEntry.DistEntryDims();
      DimVec destVec = destEntry.DistEntryDims();
      DimVec srcVec2 = srcEntry2.DistEntryDims();
      DimVec destVec2 = destEntry2.DistEntryDims();
      
      //      cout << srcEntry.PrettyStr() << " -> " << destEntry.PrettyStr() << endl;
      //      cout << srcEntry2.PrettyStr() << " -> " << destEntry2.PrettyStr() << endl;
      
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
          /*
           PrintVec(suff1);
           PrintVec(pref11);
           PrintVec(pref12);
           PrintVec(suff2);
           PrintVec(pref21);
           PrintVec(pref22);
           */
          
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
          
          //	  PrintVec(newDims1);
          intType.m_dists[m_dim].DimsToDistEntry(newDims1);
          //	  PrintVec(newDims2);
          intType.m_dists[dim].DimsToDistEntry(newDims2);
          if (DistTypeNotEqual(intType, redist->m_info.GetDist())) {
            RedistNode *newRedist = new RedistNode(intType);
            newRedist->AddInput(redist->Input(0), redist->InputConnNum(0));
            node->m_poss->AddNode(newRedist);
            
            if (intType == newRedist->InputDataType(0).GetDist())
              throw;
            
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

bool RecursivelyUpdateWithIntType(Dim minDim, Dim dim, const DistType &srcType, const DistType &fullIntType, DistType &finalIntType)
{
  if (dim >= finalIntType.m_numDims 
      || srcType.m_numDims != fullIntType.m_numDims 
      || srcType.m_numDims!= finalIntType.m_numDims) 
    {
      throw;
    }
  if (finalIntType.m_dists[dim] == fullIntType.m_dists[dim])
    return true;
      
  finalIntType.m_dists[dim] = fullIntType.m_dists[dim];

  DistEntry srcEntry = srcType.m_dists[dim];
  DistEntry newEntry = finalIntType.m_dists[dim];

  //if a grid mode disappears from src to final, then recurse on the tensor mode to which it is moved
  const DimSet srcSet = srcEntry.DistEntryDimSet();
  const DimSet newSet = newEntry.DistEntryDimSet();
  
  DimSetConstIter iter = srcSet.begin();
  for(; iter != srcSet.end(); ++iter) {
    if (newSet.find(*iter) == newSet.end()) {
      Dim dest;
      if (fullIntType.FindGridMode(*iter, dest)) {
	if (dest < minDim)
	  return false;
	if (!RecursivelyUpdateWithIntType(minDim, dest, srcType, fullIntType, finalIntType))
	  return false;
      }
    }
  }

  iter = newSet.begin();
  for(; iter != newSet.end(); ++iter) {
    if (srcSet.find(*iter) == srcSet.end()) {
      Dim dest;
      if (srcType.FindGridMode(*iter, dest)) {
	if (dest < minDim)
	  return false;
	if (!RecursivelyUpdateWithIntType(minDim, dest, srcType, fullIntType, finalIntType))
	  return false;
      }
    }
  }
  return true;
}

bool MultiIndexAllToAll::CanApply(const Node *node) const
{
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();
  
  const DimSet usedDims = destType.UsedGridDims();
  
  const Dim numDims = srcType.m_numDims;
  
  if (numDims != destType.m_numDims)
    throw;
  
  if (m_dim >= numDims)
    return false;
  
  DistType intType = srcType;
  
  DistEntryVec movedDims;
  movedDims.resize(numDims);
  DimVec disappearingDims;
  
  for (Dim dim = 0; dim < numDims; ++dim) {
    DimVec vec = intType.m_dists[dim].DistEntryDims();
    DimSet destDims = destType.m_dists[dim].DistEntryDimSet();
    //Want to find ending dims of this source mode that disappear in the
    // dest dist
    DimVecRevIter iter = vec.rbegin();
    for(; iter != vec.rend(); ++iter) {
      if (usedDims.find(*iter) != usedDims.end()) {
        break;
      }
    }
    //Now, check the last dim in this source mode that does not disappear
    //If that dim is also in the dest dim, then these disappearing dims should be handled by
    // an allgather
    //If that dim is in a different dest dim, though, then these dims should disappear
    // as part of the intermediate type so that ending dim can be part of the AllToAll
    if (iter != vec.rend()) {
      if (destDims.find(*iter) == destDims.end()) {
        //add the ending dims that disappear to disappearingDims
        // and remove them from vec so the code below can figure out the
        DimVecIter fiter = vec.begin();
        while (*fiter != *iter)
          ++fiter;
        ++fiter;
        disappearingDims.insert(disappearingDims.end(), fiter, vec.end());
        vec.erase(fiter, vec.end());
      }
    }
    else {
      continue;
    }
    
    // Now, go through dims backwards and find all dims that should be moved to other
    // modes
    while (!vec.empty()) {
      Dim currDim = vec.back();
      if (destType.m_dists[dim].ContainsDim(currDim)) {
        // this dim is still in the dest mode
        break;
      }
      else if (usedDims.find(currDim) != usedDims.end()) {
        bool found = false;
        for(Dim checkingDim = 0; !found && checkingDim < numDims; ++checkingDim) {
          if (checkingDim != dim) {
            if(destType.m_dists[checkingDim].ContainsDim(currDim)) {
              found = true;
              DistEntry tmp = movedDims[checkingDim];
              tmp.AppendDim(currDim);
              movedDims[checkingDim] = tmp;
            }
          }
        }
        if (!found)
          throw;
        vec.pop_back();
      }
      else
        break;
    }
    intType.m_dists[dim].DimsToDistEntry(vec);
  }
  
  //check modes with all dims that disappear from src, if dest has dims
  // that move to it, erase all disappearing dims from intermediate
  for(Dim dim = 0; dim < numDims; ++dim) {
    if (!movedDims[dim].IsStar()) {
      //dims are moving to this mode;
      DimVec vec = intType.m_dists[dim].DistEntryDims();
      while (!vec.empty()) {
        Dim currDim = vec.back();
        if (usedDims.find(currDim) == usedDims.end()) {
          disappearingDims.push_back(currDim);
          vec.pop_back();
        }
        else
          break;
      }
      intType.m_dists[dim].DimsToDistEntry(vec);
    }
  }
  
  //Now add moved dims to intType.
  // Add moved dims in the same order they're found in the dest
  // just in case the int dist is the same as the dest
  for(Dim dim = 0; dim < numDims; ++dim) {
    if (!movedDims[dim].IsStar()) {
      //dims are moving to this mode;
      DimVec vec = intType.m_dists[dim].DistEntryDims();
      DimSet movedSet = movedDims[dim].DistEntryDimSet();
      DimVec destVec = destType.m_dists[dim].DistEntryDims();
      DimVecIter iter = destVec.begin();
      for(; !movedSet.empty() && iter != destVec.end(); ++iter) {
        Dim currDim = *iter;
        if (movedSet.find(currDim) != movedSet.end()) {
          vec.push_back(currDim);
          movedSet.erase(currDim);
        }
      }
      if (!movedSet.empty())
        throw;
      intType.m_dists[dim].DimsToDistEntry(vec);
    }
  }

  DistType finalIntType = srcType;
  if (!RecursivelyUpdateWithIntType(m_dim, m_dim, srcType, intType, finalIntType))
    return false;

  
  return (finalIntType != destType && finalIntType != srcType);
}

void MultiIndexAllToAll::Apply(Node *node) const
{
  RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDataType(0).GetDist();
  const DistType &destType = redist->m_info.GetDist();
  
  const DimSet usedDims = destType.UsedGridDims();
  
  const Dim numDims = srcType.m_numDims;
  
  if (numDims != destType.m_numDims)
    throw;
  
  DistType intType = srcType;
  
  DistEntryVec movedDims;
  movedDims.resize(numDims);
  DimVec disappearingDims;
  
  for (Dim dim = 0; dim < numDims; ++dim) {
    DimVec vec = intType.m_dists[dim].DistEntryDims();
    DimSet destDims = destType.m_dists[dim].DistEntryDimSet();
    //Want to find ending dims of this source mode that disappear in the
    // dest dist
    DimVecRevIter iter = vec.rbegin();
    for(; iter != vec.rend(); ++iter) {
      if (usedDims.find(*iter) != usedDims.end()) {
        break;
      }
    }
    //Now, check the last dim in this source mode that does not disappear
    //If that dim is also in the dest dim, then these disappearing dims should be handled by 
    // an allgather
    //If that dim is in a different dest dim, though, then these dims should disappear
    // as part of the intermediate type so that ending dim can be part of the AllToAll
    if (iter != vec.rend()) {
      if (destDims.find(*iter) == destDims.end()) {
        //add the ending dims that disappear to disappearingDims
        // and remove them from vec so the code below can figure out the 
        DimVecIter fiter = vec.begin();
        while (*fiter != *iter)
          ++fiter;
        ++fiter;
        disappearingDims.insert(disappearingDims.end(), fiter, vec.end());
        vec.erase(fiter, vec.end());
      }
    }
    else {
      continue;
    }
    
    // Now, go through dims backwards and find all dims that should be moved to other
    // modes
    while (!vec.empty()) {
      Dim currDim = vec.back();
      if (destType.m_dists[dim].ContainsDim(currDim)) {
        // this dim is still in the dest mode
        break; 
      }
      else if (usedDims.find(currDim) != usedDims.end()) {
        bool found = false;
        for(Dim checkingDim = 0; !found && checkingDim < numDims; ++checkingDim) {
          if (checkingDim != dim) {
            if(destType.m_dists[checkingDim].ContainsDim(currDim)) {
              found = true;
              DistEntry tmp = movedDims[checkingDim];
              tmp.AppendDim(currDim);
              movedDims[checkingDim] = tmp;
            }
          }
        }
        if (!found)
          throw;
        vec.pop_back();
      }
      else
        break;
    }
    intType.m_dists[dim].DimsToDistEntry(vec);
  }
  
  //check modes with all dims that disappear from src, if dest has dims
  // that move to it, erase all disappearing dims from intermediate
  for(Dim dim = 0; dim < numDims; ++dim) {
    if (!movedDims[dim].IsStar()) {
      //dims are moving to this mode;
      DimVec vec = intType.m_dists[dim].DistEntryDims();
      while (!vec.empty()) {
        Dim currDim = vec.back();
        if (usedDims.find(currDim) == usedDims.end()) {
          disappearingDims.push_back(currDim);
          vec.pop_back();
        }
        else
          break;
      }
      intType.m_dists[dim].DimsToDistEntry(vec);
    }
  }
  
  //Now add moved dims to intType.
  // Add mvoed dims in the same order they're found in the dest
  // just in case they int dist is the same as the dest
  for(Dim dim = 0; dim < numDims; ++dim) {
    if (!movedDims[dim].IsStar()) {
      //dims are moving to this mode;
      DimVec vec = intType.m_dists[dim].DistEntryDims();
      DimSet movedSet = movedDims[dim].DistEntryDimSet();
      DimVec destVec = destType.m_dists[dim].DistEntryDims();
      DimVecIter iter = destVec.begin();
      for(; !movedSet.empty() && iter != destVec.end(); ++iter) {
        Dim currDim = *iter;
        if (movedSet.find(currDim) != movedSet.end()) {
          vec.push_back(currDim);
          movedSet.erase(currDim);
        }
      }
      if (!movedSet.empty())
        throw;
      intType.m_dists[dim].DimsToDistEntry(vec);
    }
  }

  DistType finalIntType = srcType;
  if (!RecursivelyUpdateWithIntType(m_dim, m_dim,srcType, intType, finalIntType))
    throw;

  RedistNode *intRedist = new RedistNode(finalIntType);
  intRedist->AddInput(redist->Input(0), redist->InputConnNum(0));
  
  RedistNode *finalRedist = new RedistNode(destType, redist->m_info.GetPerm(),
                                           redist->m_align, redist->m_alignModes, redist->m_alignModesSrc);
  finalRedist->AddInput(intRedist, 0);
  
  redist->RedirectAllChildren(finalRedist);
  
  redist->m_poss->AddNode(intRedist);
  redist->m_poss->AddNode(finalRedist);
  
  redist->m_poss->DeleteChildAndCleanUp(redist);
}

bool FindPermCycle(const DistType &srcType, const DistType &destType,
		   DimSet &tensorModesInvolved, 
		   const DimVec &gridModesToFind,
		   Size suffixLen,
		   DimVec *gridModesInvolved)
{
  for (Dim dim = 0; dim < srcType.m_numDims; ++dim) {
    DistEntry srcDistEntry = srcType.m_dists[dim];
    DistEntry destDistEntry = destType.m_dists[dim];
    if (srcDistEntry != destDistEntry) {
      DimVec srcDims = srcDistEntry.DistEntryDims();
      DimVec destDims = destDistEntry.DistEntryDims();
      DimVec pref;
      DimVec srcSuff, destSuff;
      GetCommonPrefix(srcDims, destDims,
		      pref,
		      srcSuff, destSuff);
      if (!destSuff.empty()
          && std::is_permutation(destSuff.begin(), destSuff.end(),
			      gridModesToFind.begin())) 
	{
	  if (tensorModesInvolved.find(dim) != tensorModesInvolved.end()) {
	    return true;
	  }
	  //If we're not allowing the optimization of use permutation
	  // for moving (permuted) suffixes around if the suffixes have
	  // the same number of processes, then recursion isn't allowed
	  // since we only allow permutation on the single tensor mode
	  if (!M_allowSquareGridOpt)
	    return false;
	  tensorModesInvolved.insert(dim);
	  if (gridModesInvolved) {
	    gridModesInvolved->insert(gridModesInvolved->end(),
				      srcSuff.begin(), srcSuff.end());
	  }
    if ((srcSuff.empty() ? 0 : GridModeLens(srcSuff)) != suffixLen)
	    return false;
	  return FindPermCycle(srcType, destType,
			       tensorModesInvolved,
			       srcSuff,
			       suffixLen,
			       gridModesInvolved);
	}
    }
  }
  return false;
}

bool GetPermPattern(const DistType &srcType,
                         const DistType &destType,
                         DimVec *gridModes)
{
  Dim numDims = srcType.m_numDims;
  
  if (numDims != destType.m_numDims)
    return false;
  
  if (gridModes)
    gridModes->clear();

  DimSet tensorModesInvolved;
  
  for (Dim dim = 0; dim < numDims; ++dim) {
    if (tensorModesInvolved.find(dim) == tensorModesInvolved.end())
      {
	DistEntry srcDistEntry = srcType.m_dists[dim];
	DistEntry destDistEntry = destType.m_dists[dim];
	if (srcDistEntry != destDistEntry) {
	  DimVec srcDims = srcDistEntry.DistEntryDims();
	  DimVec destDims = destDistEntry.DistEntryDims();
	  DimVec pref;
	  DimVec srcSuff, destSuff;
	  GetCommonPrefix(srcDims, destDims,
			  pref,
			  srcSuff, destSuff);
	  
	  if (!srcSuff.empty() && !destSuff.empty()) {
	    Size suffixLen = GridModeLens(srcSuff);
	    if (suffixLen != GridModeLens(destSuff))
	      return false;
	    DimSet localTensorModesInvolved;
	    localTensorModesInvolved.insert(dim);
	    if (gridModes) {
	      gridModes->insert(gridModes->end(),
					srcSuff.begin(), srcSuff.end());
	    }	    
	    if (!FindPermCycle(srcType, destType,
			       localTensorModesInvolved,
			       srcSuff,
			       suffixLen,
			       gridModes))
	      {
		if (gridModes)
		  gridModes->clear();
		return false;
	      }
	    else {
	      tensorModesInvolved.insert(localTensorModesInvolved.begin(), localTensorModesInvolved.end());
	    }
	  }
	  else
	    return false;
	}
      }
  }

  return !tensorModesInvolved.empty();
}

string GetAlignmentSource(Node *node, ConnNum inNum)
{
  Node *in = node->Input(inNum);
  if (in->GetNodeClass() == TempVarNode::GetClass()) {
    TempVarNode *tmp = (TempVarNode*)in;
    if (tmp->m_sumDims.empty())
      return GetAlignmentSource(tmp, 0);
  }
  return in->GetNameStr(node->InputConnNum(inNum));
}
#endif



