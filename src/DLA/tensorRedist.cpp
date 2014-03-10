#include "tensorRedist.h"

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
  if (m_inputs.empty()) {
    return "RedistNode to " + DistTypeToStr(m_destType) +
      " without parent";
  }
  else {
    Node *parent = Input(0);
    DistType type = ((DLANode*)parent)->GetDistType(InputConnNum(0));
    return  "RedistNode " + DistTypeToStr(type) +
      " -> " + DistTypeToStr(m_destType);
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
    DLANode *parent = (DLANode*)Input(0);
    DistType m_srcType = parent->GetDistType(InputConnNum(0));
    parent->Prop();
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

#endif
