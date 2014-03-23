#include "tensorRedist.h"
#include <ostream>
#include <sstream>

#if DOTENSORS

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
    delete [] m_lsizes;
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
    
    if (numDims != InputNumDims(0))
      throw;

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
      m_cost = 0;
    }
    else if (diffs.size() == 1) {
      const Dim dim = *(diffs.begin());
      DimVec src = DistType::DistEntryDims(m_srcType.m_dists[dim]);
      DimVec dest = DistType::DistEntryDims(m_destType.m_dists[dim]);

      
      if (src.empty() || IsPrefix(src, dest)) {
	//local memory copy
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
    else {
      m_cost = 0;
      const unsigned int totNumIters = m_lsizes[0].NumSizes();
      unsigned int numProcs = 1;
      DimSet unionSet;

      DimSetIter diffIter = diffs.begin();
      for(; diffIter != diffs.end(); ++diffIter) {
	Dim diffDim = *diffIter;
	DimVec src = DistType::DistEntryDims(m_srcType.m_dists[diffDim]);
	DimVec dest = DistType::DistEntryDims(m_destType.m_dists[diffDim]);

	if (src.empty() || IsPrefix(src, dest)) {
	  //local mem copy for this dimensions
	}
	else if (IsPrefix(dest, src) || dest.empty()) {
	  DimVecIter iter = src.begin() + dest.size();
	  for(; iter != src.end(); ++iter) {
	    if (unionSet.insert(*iter).second)
	      numProcs *= GridLens[*iter];
	  }
	}
	else {
	  DimVecIter iter = src.begin();
	  for(; iter != src.end(); ++iter) {
	    if (unionSet.insert(*iter).second)
	      numProcs *= GridLens[*iter];
	  }
	  iter = dest.begin();
	  for(; iter != dest.end(); ++iter) {
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
    throw;
  }
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
  return m_lsizes+dim;
}

void RedistNode::ClearSizeCache()
{
  if (m_lsizes) {
    delete m_lsizes;
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
  m_lsizes = new Sizes[numDims];
  for (Dim dim = 0; dim < numDims; ++dim)
    GetLocalSizes(m_destType, dim, in->Len(num,dim), m_lsizes+dim);
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
  out.Indent();
  DistType in = InputDistType(0);

  *out << GetName(0).str()
       << " = "
       << GetInputName(0).str()
       << ";\n";
}


RedistNodeWithSummation::RedistNodeWithSummation(const DistType &destType, const DimVec &sumDims)
  : RedistNode(destType)
{
  m_sumDims = sumDims;
}


void RedistNodeWithSummation::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const RedistNodeWithSummation *node = (RedistNodeWithSummation*)orig;
  RedistNode::Duplicate(node, shallow, possMerging);
  m_sumDims = node->m_sumDims;  
}

NodeType RedistNodeWithSummation::GetType() const
{
  stringstream str;
  str << "redistWithSum";
  DimVecConstIter iter = m_sumDims.begin();
  for(; iter != m_sumDims.end(); ++iter)
    str << *iter << ",";
  return str.str();
}
void RedistNodeWithSummation::SanityCheck()
{
  throw;
}

void RedistNodeWithSummation::Prop()
{
  throw;
}

void RedistNodeWithSummation::PrintCode(IndStream &out)
{
  throw;
}

void RedistNodeWithSummation::FlattenCore(ofstream &out) const
{
  throw;
}

void RedistNodeWithSummation::UnflattenCore(ifstream &in, SaveInfo &info)
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
#endif
