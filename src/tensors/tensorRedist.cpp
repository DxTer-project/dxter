#include "tensorRedist.h"
#include <ostream>
#include <sstream>

#if DOTENSORS

#include "tensorSumScatter.h"

void GetCommonSuffix(const DimVec &dims1, const DimVec &dims2, 
		     DimVec &suff,
		     DimVec &pref1, DimVec &pref2);
void GetDifferentSuffix(const DimVec &dims1, const DimVec &dims2, 
			DimVec &suff1, DimVec &suff2);
void GetSuffix(const DimVec &dims1, const DimVec &dims2, 
	       DimVec &suff);

RedistNode::RedistNode() 
{
  m_destType.SetToDefault(0);
  m_lsizes = NULL;
}

RedistNode::RedistNode(const DistType &destType)
{
  m_destType = destType;
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
  m_destType = origNode->m_destType;
}

NodeType RedistNode::GetType() const
{
  if (m_name.length()) 
    return m_name;
  else {
    return (string)"RedistNode to " +  m_destType.QuickStr();
  }
}

void RedistNode::SanityCheck()
{
  DLANode::SanityCheck();
  if (m_inputs.size() != 1) {
    cout << "m_inputs.size() != 1\n";
    throw;
  }
  
  if (!m_children.size())
    throw;
}

void RedistNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (!m_name.length())
      m_name = (string)"RedistNode to " +  m_destType.QuickStr();
    DLANode *parent = (DLANode*)Input(0);
    const DistType &m_srcType = parent->GetDistType(InputConnNum(0));
    parent->Prop();
    const Dim numDims = m_destType.m_numDims;
    
    if (numDims != InputNumDims(0)) {
      cout << "numDims " << numDims << " " << InputNumDims(0) << endl;
      cout << "input " << Input(0)->GetNodeClass() << endl;
      cout << GetInputNameStr(0) << endl;
      cout << m_destType.PrettyStr() << endl;
      throw;
    }
    if (m_srcType.m_numDims != numDims)
      throw;

    if (m_srcType.m_numDims != numDims) {
      cout << m_destType.str() << " <- " << m_srcType.str() << endl;
      throw;
    }

    if (!m_destType.IsSane()) {
      cout << m_destType.str() << endl;
      m_poss->PrintTransVec();
      throw;
    }

    DimSet diffs;
    
    for (Dim dim = 0; dim < numDims; ++dim) {
      if (m_srcType.m_dists[dim] != m_destType.m_dists[dim]) {
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
      DimVec dest = m_destType.m_dists[dim].DistEntryDims();

      
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
      m_cost = 0;
      const unsigned int totNumIters = m_lsizes[0].NumSizes();
      unsigned int numProcs = 1;
      DimSet unionSet;
      
      DimSetIter diffIter = diffs.begin();
      for(; diffIter != diffs.end(); ++diffIter) {
	Dim diffDim = *diffIter;
	DimVec src = m_srcType.m_dists[diffDim].DistEntryDims();
	DimVec dest = m_destType.m_dists[diffDim].DistEntryDims();

	if (src.empty() || IsPrefix(src, dest)) {
	  m_cost = 0;
	  //local mem copy for this dimensions
	}
	else if (IsPrefix(dest, src) || dest.empty()) {
	  m_cost = 0;
	  //	  cout << "AllGather with AllToAll\n";
	  //	  throw;
	  //	  DimVecIter iter = src.begin() + dest.size();
	  //	  for(; iter != src.end(); ++iter) {
	  //	    if (unionSet.insert(*iter).second)
	  //	      numProcs *= GridLens[*iter];
	  //	  }
	}
	else {
	  DimVec suff1, suff2;
	  GetDifferentSuffix(src, dest, suff1, suff2);
	  DimVecIter iter = suff1.begin();
	  for(; iter != suff1.end(); ++iter) {
	    if (unionSet.insert(*iter).second)
	      numProcs *= GridLens[*iter];
	  }
	  iter = suff2.begin();
	  for(; iter != suff2.end(); ++iter) {
	    if (unionSet.insert(*iter).second)
	      numProcs *= GridLens[*iter];
	  }
	}
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
  const DistType &m_srcType = parent->GetDistType(InputConnNum(0));
  const Dim numDims = m_destType.m_numDims;
  
  bool foundMemCopy = false;
  bool foundAllGather = false;
  bool foundAllToAll = false;
  DimSet diffs;
  
  for (Dim dim = 0; dim < numDims; ++dim) {
    DistEntry srcDistEntry = m_srcType.m_dists[dim];
    DistEntry destDistEntry = m_destType.m_dists[dim];
    if (srcDistEntry != destDistEntry) {
      DimVec srcDims = srcDistEntry.DistEntryDims();
      DimVec destDims = destDistEntry.DistEntryDims();
      if (srcDistEntry.IsStar() || IsPrefix(srcDims, destDims)) {
	foundMemCopy = true;
	if (foundAllGather || foundAllToAll) {
	  return ROTENSORPHASE;
	}
      }
      else if (destDistEntry.IsStar() || IsPrefix(destDims, srcDims)) {
	if (foundAllGather || foundMemCopy || foundAllToAll)
	  return ROTENSORPHASE;
	else
	  foundAllGather = true;
      }
      else {
	if (foundMemCopy || foundAllGather)
	  return ROTENSORPHASE;
	else {
	  foundAllToAll = true;
	  if (diffs.size() > 2)
	    return ROTENSORPHASE;
	  else if (diffs.size() == 1) {
	    Dim otherDim = *(diffs.begin());
	    DimVec otherSrcDims = m_srcType.m_dists[otherDim].DistEntryDims();
	    DimVec otherDestDims = m_destType.m_dists[otherDim].DistEntryDims();
	    
	    DimVec commonSuff1, commonSuff2, otherSrcPref, destPref, srcPref, otherDestPref;
	    GetCommonSuffix(otherSrcDims, destDims, commonSuff1, otherSrcPref, destPref);
	    GetCommonSuffix(srcDims, otherDestDims, commonSuff2, srcPref, otherDestPref);
	    
	    if (commonSuff1.empty())
	      return ROTENSORPHASE;

	    if (commonSuff1 != commonSuff2) {
	      return ROTENSORPHASE;
	    }
	    if (otherSrcPref != otherDestPref) {
	      return ROTENSORPHASE;
	    }
	    if (srcPref != destPref) {
	      return ROTENSORPHASE;
	    }
	  }
	  diffs.insert(dim);
	}
      }
    }
  }

  if (diffs.size() == 1 && foundAllToAll) {
    Dim dim = *(diffs.begin());
    DistEntry srcDistEntry = m_srcType.m_dists[dim];
    DistEntry destDistEntry = m_destType.m_dists[dim];

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

void RedistNode::ClearSizeCache()
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

void RedistNode::BuildSizeCache()
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
      GetLocalSizes(m_destType, dim, in->Len(num,dim), m_lsizes+dim);
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
  name.m_type = m_destType;
  return name;
}

void RedistNode::PrintCode(IndStream &out)
{  
  //Reflect in AddVars
  out.Indent();

  const DistType &m_srcType = InputDistType(0);
  const Dim numDims = m_destType.m_numDims;

  string inName = GetInputName(0).str();
  string outName = GetName(0).str(); 

  *out << "   // " << GetName(0).PrettyStr() 
       << " <- " << GetInputName(0).PrettyStr() << endl;
  out.Indent();
  
  if (m_srcType.m_numDims != numDims)
    throw;

  DimSet diffs;

  *out << outName << ".SizeTo( " << inName << " );\n";
  out.Indent();
    
  for (Dim dim = 0; dim < numDims; ++dim) {
    if (m_srcType.m_dists[dim] != m_destType.m_dists[dim]) {
      diffs.insert(dim);
    }
  }

  if (diffs.empty()) {
    throw;
  }
  else if (diffs.size() == 1) {
    const Dim dim = *(diffs.begin());
    DimVec src = m_srcType.m_dists[dim].DistEntryDims();
    DimVec dest = m_destType.m_dists[dim].DistEntryDims();
       
    if (src.empty() || IsPrefix(src, dest)) {
      DimVec suff;
      GetSuffix(src, dest, suff);
      *out << "LocalRedist( " << outName << ", "
	   << inName << ", " << dim 
	   << ", " << ModeArrayVarName(suff)
	   << " );\n";
    }
    else if (IsPrefix(dest, src) || dest.empty()) {
      if (!src.empty()) {
	
	*out << "AllGatherRedistPartial( " << outName << ", "
	     << inName << ", " << dim << ", ";
	DimVec suff;
	GetSuffix(dest, src, suff);
	*out << ModeArrayVarName(suff) << " );\n";

      }
      else {
	*out << "AllGatherRedist( " << outName << ", "
	     << inName << ", " << dim << " );\n";
      }
    }
    else {
      *out << "PermutationRedist( " << outName << ", "
      	   << inName << ", " << dim << " );\n";
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

    *out << "AllToAllDoubleIndexRedist( "  << outName << ", "
	 << inName << ", " 
	 << IndexPairVarName(diff1, diff2)
	 << ", ";

    
    DimVec firstVec, secondVec;
    bool first = true;
    iter = diffs.begin();
    for(; iter != diffs.end(); ++iter) {
      Dim diffDim = *iter;
      DimVec src = m_srcType.m_dists[diffDim].DistEntryDims();
      DimVec dest = m_destType.m_dists[diffDim].DistEntryDims();

      if (src.empty() || IsPrefix(src, dest)) {
	throw;
      }
      else if (IsPrefix(dest, src) || dest.empty()) {
	throw;
      }
      else {
	DimVec suff1, suff2;
	GetDifferentSuffix(src, dest, suff1, suff2);
	if (first)
	  firstVec = suff2;
	else
	  secondVec = suff2;
      }
      if (first) {
	first = false;
      }
    }

    *out << ModeArrayPairVarName(firstVec, secondVec)
	 << " );\n";
  }
}

void RedistNode::AddVariables(VarSet &set) const
{
  //Reflect in PrintCode
  DLANode::AddVariables(set);
  
  const DistType &m_srcType = InputDistType(0);
  const Dim numDims = m_destType.m_numDims;

  string inName = GetInputName(0).str();
  string outName = GetName(0).str(); 

  if (m_srcType.m_numDims != numDims)
    throw;

  DimSet diffs;
    
  for (Dim dim = 0; dim < numDims; ++dim) {
    if (m_srcType.m_dists[dim] != m_destType.m_dists[dim]) {
      diffs.insert(dim);
    }
  }

  if (diffs.empty()) {
    throw;
  }
  else if (diffs.size() == 1) {
    const Dim dim = *(diffs.begin());
    DimVec src = m_srcType.m_dists[dim].DistEntryDims();
    DimVec dest = m_destType.m_dists[dim].DistEntryDims();
       
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
    bool first = true;
    iter = diffs.begin();
    for(; iter != diffs.end(); ++iter) {
      Dim diffDim = *iter;
      DimVec src = m_srcType.m_dists[diffDim].DistEntryDims();
      DimVec dest = m_destType.m_dists[diffDim].DistEntryDims();

      if (src.empty() || IsPrefix(src, dest)) {
	throw;
      }
      else if (IsPrefix(dest, src) || dest.empty()) {
	throw;
      }
      else {
	DimVec suff1, suff2;
	GetDifferentSuffix(src, dest, suff1, suff2);
	if (first)
	  firstVec = suff2;
	else
	  secondVec = suff2;
      }
      if (first) {
	first = false;
      }
    }
    Var arrayPairVar(firstVec, secondVec);
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
void AllReduceNode::SanityCheck()
{
  throw;
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



bool RemoveWastedRedist::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  RedistNode *redistNode = (RedistNode*)node;
  const DistType* type = &(redistNode->m_destType);
  if (node->m_children.size() == 0)
    throw;
  while(redistNode->Input(0) 
        && (redistNode->Input(0)->GetNodeClass() == RedistNode::GetClass()))
    {
      redistNode = (RedistNode*)redistNode->Input(0);
      if (DistTypeEqual(redistNode->m_destType,*type))
	return true;
      for(unsigned int i = 0; i < redistNode->m_children.size(); ++i) {
	Node *tmp = redistNode->Child(i);
	if (tmp != node && tmp->GetNodeClass() == RedistNode::GetClass()) {
	  if (DistTypeEqual(((RedistNode*)tmp)->m_destType, *type))
	    return true;
	}
      }
    }
  if (redistNode->Input(0)
      && DistTypeEqual(((DLANode*)(redistNode->Input(0)))->GetDistType(redistNode->InputConnNum(0)),*type))
    return true;
  return false;
}

void RemoveWastedRedist::Apply(Poss *poss, Node *node) const
{
  RedistNode *redistNode = (RedistNode*)node;
  const DistType* type = &(redistNode->m_destType);
  while(redistNode->Input(0) 
        && (redistNode->Input(0)->GetNodeClass() == RedistNode::GetClass()))
    {
      redistNode = (RedistNode*)redistNode->Input(0);
      if (DistTypeEqual(redistNode->m_destType, *type)) {
	node->RedirectChildren(redistNode, 0);
	node->m_poss->DeleteChildAndCleanUp(node);
	return;
      }
      for(unsigned int i = 0; i < redistNode->m_children.size(); ++i) {
	Node *tmp = redistNode->Child(i);
	if (tmp != node && tmp->GetNodeClass() == RedistNode::GetClass()) {
	  if (DistTypeEqual(((RedistNode*)tmp)->m_destType,*type)) {
	    node->RedirectChildren(tmp, 0);
	    node->m_poss->DeleteChildAndCleanUp(node);
	    return;
	  }
	}
      }
    }
  if (redistNode->Input(0)
      && DistTypeEqual(((DLANode*)redistNode->Input(0))->GetDistType(redistNode->InputConnNum(0)),*type))
    {
      node->RedirectChildren(redistNode->Input(0), redistNode->InputConnNum(0));
      node->m_poss->DeleteChildAndCleanUp(node);
      return;
    }
  throw;
}


bool RemoveNOPRedistribs::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const DLANode *ddla = (DLANode*)node;
  const Node *parent = ddla->Input(0);
  if (!parent
      || (DistTypeNotEqual(((DLANode*)parent)->GetDistType(node->InputConnNum(0)), ddla->GetDistType(0)))) {
    return false;
  }
  return true;
}

void RemoveNOPRedistribs::Apply(Poss *poss, Node *node) const
{
  Node *parent = node->Input(0);
  node->RedirectChildren(parent,node->InputConnNum(0));
  node->m_poss->DeleteChildAndCleanUp(node);
}

bool CombineRedistribs::CanApply(const Poss *poss, const Node *node) const
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
        && DistTypeEqual(((RedistNode*)output)->m_destType, redist->m_destType)) 
      {
	return true;
      }
  }
  return false;
}

void CombineRedistribs::Apply(Poss *poss, Node *node) const
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
        && DistTypeEqual(((RedistNode*)output)->m_destType, redist->m_destType)) 
      {
	output->RedirectChildren(node,0);
	output->m_poss->DeleteChildAndCleanUp(output);
	return;
      }
  }
  throw;
}

bool SplitRedistribs::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &src = redist->InputDistType(0);
  const DistType *dest = &(redist->m_destType);

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

void SplitRedistribs::Apply(Poss *poss, Node *node) const
{
  RedistNode *orig = (RedistNode*)node;
  DistType one = orig->InputDistType(0);
  const DistType *two = &(orig->m_destType);

  one.m_dists[m_dim] = two->m_dists[m_dim];

/*  if (!one.IsSane()) {
    cout << "splitting " << m_dim << endl;
    cout << DistTypeToStr(one) << endl;
    cout << "from " << DistTypeToStr(orig->InputDistType(0))
	 << " -> " << DistTypeToStr(orig->m_destType);
    throw;
    }*/



  RedistNode *newRedist = new RedistNode(one);
  newRedist->AddInput(orig->Input(0), orig->InputConnNum(0));
  poss->AddNode(newRedist);

  RedistNode *newRedist2 = new RedistNode(*two);
  newRedist2->AddInput(newRedist, 0);
  poss->AddNode(newRedist2);

  node->RedirectChildren(newRedist2, 0);
  node->m_poss->DeleteChildAndCleanUp(node);

}

bool SingleIndexAllToAll::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDistType(0);

  if (srcType.m_numDims != redist->m_destType.m_numDims)
    throw;
  else if (srcType.m_numDims <= m_dim)
    return false;

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = redist->m_destType.m_dists[m_dim];
  if (srcEntry == destEntry)
    return false;

  for(Dim dim = 0; dim < srcType.m_numDims; ++dim) {
    if (dim != m_dim) {
      if (srcType.m_dists[dim] != redist->m_destType.m_dists[dim])
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

void SingleIndexAllToAll::Apply(Poss *poss, Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDistType(0);

  DistType type1 = srcType;
  DistEntry entry1 = type1.m_dists[m_dim];
  DimVec vec1 = entry1.DistEntryDims();
  DimSet set1;
  set1.insert(vec1.begin(), vec1.end());
  
  DimVec destVec = redist->m_destType.m_dists[m_dim].DistEntryDims();
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
  poss->AddNode(redist1);

  DistType type2 = redist->m_destType;
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
  poss->AddNode(redist2);

  RedistNode *redist3 = new RedistNode(redist->m_destType);
  redist3->AddInput(redist2, 0);
  poss->AddNode(redist3);

  node->RedirectChildren(redist3, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

bool SplitAllGathers::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDistType(0);
  const DistType &destType = redist->m_destType;

  if (srcType.m_numDims != redist->m_destType.m_numDims)
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

void SplitAllGathers::Apply(Poss *poss, Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  const RedistNode *redist = (RedistNode*)node;
  const DistType &srcType = redist->InputDistType(0);
  const DistType &destType = redist->m_destType;

  DistEntry srcEntry = srcType.m_dists[m_dim];
  DistEntry destEntry = destType.m_dists[m_dim];

  DimVec destDims = destEntry.DistEntryDims();
  DimVec srcDims = srcEntry.DistEntryDims();

  srcDims.pop_back();

  DistType intType = srcType;
  intType.m_dists[m_dim].DimsToDistEntry(srcDims);

  RedistNode *redist1 = new RedistNode(intType);
  redist1->AddInput(redist->Input(0), redist->InputConnNum(0));
  poss->AddNode(redist1);

  RedistNode *redist2 = new RedistNode(destType);
  redist2->AddInput(redist1, 0);
  poss->AddNode(redist2);

  node->RedirectChildren(redist2, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}


#endif



