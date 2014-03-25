/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2013, The University of Texas and Bryan Marker

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



#include "helperNodes.h"
#include "elemRedist.h"
#include <cmath>


InputNode::InputNode() 
  : 
  m_type("InputNode"), 
#if TWOD
    m_mlsize(NULL), m_nlsize(NULL) 
#else
  m_numDims(0)
#endif
    {}

#if TWOD
InputNode::InputNode(NodeType type, Size m, Size n, string name)
: 
#if DODM
m_type(type),
#endif
m_msize(NAN), m_nsize(NAN), m_mlsize(NULL), m_nlsize(NULL)
{
  m_msize.AddRepeatedSizes(m, 1, 1);
  m_nsize.AddRepeatedSizes(n, 1, 1);
  m_varName.m_name = name;
#if DOELEM
  m_varName.m_type = D_MC_MR;
#endif
}
#endif

#if DODM&&TWOD
InputNode::InputNode(NodeType type, Size m, Size n, string name, DistType dist)
: m_type(type),
 m_msize(NAN), m_nsize(NAN), m_mlsize(NULL), m_nlsize(NULL)
{
  m_msize.AddRepeatedSizes(m,1,1);
  m_nsize.AddRepeatedSizes(n,1,1);
  m_varName.m_name = name;
  m_varName.m_type = dist;
}
#endif


#if DOTENSORS
InputNode::InputNode(NodeType type, const SizesArray sizes, string name, string indices)
: 
  m_type(type), m_numDims(indices.length()), m_lsizes(NULL)
{
  if (m_numDims > NUM_GRID_DIMS)
    throw;
  m_sizes = new Sizes[m_numDims];
  for(unsigned int i = 0; i < m_numDims; ++i)
    m_sizes[i] = sizes[i];
  m_varName.m_name = name;
  m_varName.m_indices = indices;
  m_varName.m_type.SetToDefault(m_numDims);
}

InputNode::InputNode(NodeType type, const SizesArray sizes, const DistType &dist, string name, string indices)
:
  m_type(type), m_numDims(dist.m_numDims), m_lsizes(NULL)
{
  m_sizes = new Sizes[dist.m_numDims];
  for(unsigned int i = 0; i < m_numDims; ++i)
    m_sizes[i] = sizes[i];
  m_varName.m_name = name;
  m_varName.m_type = dist;
  m_varName.m_indices = indices;
  if (m_varName.m_indices.size() != m_numDims)
    throw;
}
#endif

InputNode::~InputNode()
{
#if TWOD
  if (m_mlsize) {
    delete m_mlsize;
    delete m_nlsize;
  }
#else
  delete [] m_sizes;
  if (m_lsizes)
    delete [] m_lsizes;
#endif
}

void InputNode::PrintCode(IndStream &out)
{
#if DOTENSORS
  out.Indent();
  *out << "// " << m_type << " has " << m_numDims 
       << " dims and indices " << m_varName.m_indices << endl;
  out.Indent();
  *out << "//\tStarting distribution: " << DistTypeToStr(m_varName.m_type) << endl;
#endif
}

void InputNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const InputNode *node = (InputNode*)orig;
  m_type = node->m_type;
#if TWOD
  m_msize = node->m_msize;
  m_nsize = node->m_nsize;
#else
  m_numDims = node->m_numDims;
  m_sizes = new Sizes[m_numDims];
  for (unsigned int i = 0; i < m_numDims; ++i)
    m_sizes[i] = node->m_sizes[i];
  if(node->m_lsizes) {
    m_lsizes = new Sizes[m_numDims];
    for (unsigned int i = 0; i < m_numDims; ++i)
      m_lsizes[i] = node->m_lsizes[i];
  }
#endif
  m_varName = node->m_varName;
}

void InputNode::SanityCheck()
{
  DLANode::SanityCheck();
  if (!m_inputs.empty()) {
    cout << "!m_inputs.empty()\n";
    throw;
  }
}

void InputNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    m_cost = ZERO;
  }
}

#if TWOD
const Sizes* InputNode::GetM(unsigned int num) const
{
  if (num > 0)
    throw;
  return &m_msize;
}

const Sizes* InputNode::GetN(unsigned int num) const
{
  if (num > 0)
    throw;
  return &m_nsize;
}


const Sizes* InputNode::LocalM(unsigned int num) const
{
  if (num > 0)
    throw;
  return m_mlsize;
}

const Sizes* InputNode::LocalN(unsigned int num) const
{
  if (num > 0)
    throw;
  return m_nlsize;
}
#else

const Dim InputNode::NumDims(unsigned int num) const
{
  if (num > 0)
    throw;
  return m_numDims;
}


const Sizes* InputNode::Len(unsigned int num,Dim dim) const
{
  if (num > 0)
    throw;
  if (dim >= m_numDims)
    throw;
  return m_sizes+dim;
}


const Sizes* InputNode::LocalLen(unsigned int num,Dim dim) const
{
  if (num > 0)
    throw;
  if (dim >= m_numDims)
    throw;
  if (!m_lsizes)
    throw;
  return m_lsizes+dim;
}
#endif

Name InputNode::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  return m_varName;
}

#if TWOD

void InputNode::ClearSizeCache()
{
  if (!m_mlsize)
    return;
  delete m_mlsize;
  m_mlsize = NULL;
  delete m_nlsize;
  m_nlsize = NULL;
}
#else

void InputNode::ClearSizeCache()
{
  if (m_lsizes)
    delete [] m_lsizes;
  m_lsizes = NULL;
}
#endif

#if TWOD
void InputNode::BuildSizeCache()
{
  if (m_mlsize)
    return;
  m_mlsize = new Sizes;
  m_nlsize = new Sizes;
#if DODM
  GetLocalSizes(m_varName.m_type, &m_msize, &m_nsize, *m_mlsize, *m_nlsize);
#elif DOBLIS
  *m_mlsize = m_msize;
  *m_nlsize = m_nsize;
#else
  lkjsdf
#endif
}
#else
void InputNode::BuildSizeCache()
{
  if (m_lsizes)
    return;
  m_lsizes = new Sizes[m_numDims];
  GetLocalSizes(m_varName.m_type, m_sizes, m_lsizes);
}
#endif

void InputNode::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  WRITE(m_type);
#if TWOD
  Size size = m_msize[0];
  WRITE(size);
  size = m_nsize[0];
  WRITE(size);
  m_varName.Flatten(out);
#else
  throw;
#endif
}


void InputNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLANode::UnflattenCore(in, info);
#if TWOD
  READ(m_type);
  Size size;
  READ(size);
  m_msize.AddRepeatedSizes(size,1,1);
  READ(size);
  m_nsize.AddRepeatedSizes(size,1,1);
  m_varName.Unflatten(in);
#else
  throw;
#endif
}

void OutputNode::Duplicate(const Node *orig,bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  m_type = ((OutputNode*)orig)->m_type;
}

#if DODM
const DistType& OutputNode::GetDistType(unsigned int num) const 
{ 
#if DOELEM
  return MC_MR; 
#elif DOTENSORS
  throw;
  //  return DEFAULTDISTTYPE;
#else
  throw;
#endif
}
#endif

void OutputNode::SanityCheck()
{
  DLANode::SanityCheck();
  if (m_inputs.size() != 1) {
    cout << "m_inputs.size() != 1\n";
    throw;
  }
}

void OutputNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    Node *input = Input(0);
    input->Prop();
    m_cost = ZERO;
  }
}

#if TWOD
const Sizes* OutputNode::GetM(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputM(0);
}

const Sizes* OutputNode::GetN(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputN(0);
}


const Sizes* OutputNode::LocalM(unsigned int num) const
{
  if (num > 0)
    throw;
  return InputLocalM(0);
}

const Sizes* OutputNode::LocalN(unsigned int num) const
{
  if (num > 0)
    throw;
  return InputLocalN(0);
}

#else

const Dim OutputNode::NumDims(unsigned int num) const
{
  if (num > 0)
    throw;
  return InputNumDims(0);
}
const Sizes* OutputNode::Len(unsigned int num, Dim dim) const
{
  if (num > 0)
    throw;
  return InputLen(0,dim);
}


const Sizes* OutputNode::LocalLen(unsigned int num, Dim dim) const
{
  if (num > 0)
    throw;
  return InputLocalLen(0, dim);
}
#endif

Name OutputNode::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputName(0);
}

void OutputNode::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  out << m_type << endl;
}

void OutputNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLANode::UnflattenCore(in, info);
  getline(in, m_type);
}

#if TWOD

ConstVal::ConstVal(string name, Coef val)
  :m_val(val)
{
  m_val = val;
  m_varName.m_name = name;
#if DOELEM
  m_varName.m_type = D_STAR_STAR;
#elif DOTENSORS
  m_varName.m_type = DEFAULTDISTTYPE;
#endif
}

void ConstVal::PrintCode(IndStream &out)
{
  out.Indent();
  *out << GetName(0).str() << ".LocalBuffer(0,0) = ";
  out << m_val;
  *out << endl;
}

void ConstVal::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const ConstVal *node = (ConstVal*)orig;
  m_val = node->m_val;
  m_varName = node->m_varName;
}

void ConstVal::SanityCheck()
{
  DLANode::SanityCheck();
  if (!m_inputs.empty()) {
    cout << "!m_inputs.empty()\n";
    throw;
  }
}

void ConstVal::Prop()
{
  if (!IsValidCost(m_cost)) {
    m_cost = ZERO;
  }
}

Name ConstVal::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  return m_varName;
}

void ConstVal::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  m_varName.Flatten(out);
  WRITE(m_val);
}


void ConstVal::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  DLANode::UnflattenCore(in, info);
  m_varName.Unflatten(in);
  READ(m_val);
}
#endif

NodeType TempVarNode::GetType() const 
{
  if (m_inputs.size() != 1) {
    cout << "m_inputs.size() != 1\n";
    cout.flush();
    throw;
  }
  return GetInputNameStr(0) + " temp";
}

void TempVarNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig,shallow, possMerging);
#if DODM
  m_distType = ((TempVarNode*)orig)->m_distType;
#endif
  m_name = ((TempVarNode*)orig)->m_name;
}

void TempVarNode::PrintCode(IndStream &out)
{
  /*
  if (GetLayer() == SQ2LAYER || GetLayer() == SQ1LAYER) {
    string name = GetInputNameStr(0);
    out.Indent();
    *out << "bli_obj_create( BLIS_DOUBLE, bli_obj_length("
	 << name << "), bli_obj_width(" << name << "), 0, 0, &"
	 << GetNameStr(0) << " );\n";
    out.Indent();
    *out << "bli_copym(&" << name << ", &" << GetNameStr(0) << ");\n";
  }*/
}

void TempVarNode::SanityCheck()
{
  DLANode::SanityCheck();
  if (m_inputs.size() != 1) {
    cout << "m_inputs.size() != 1\n";
    throw;
  }
}

void TempVarNode::Prop()
{
  if (m_inputs.size() != 1) {
    cout << "m_inputs.size() != 1\n";
    cout.flush();
    throw;
  }
  if (!IsValidCost(m_cost)) {
    Input(0)->Prop();
    m_cost = ZERO;
  }
}

#if TWOD
const Sizes* TempVarNode::GetM(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputM(0);
}

const Sizes* TempVarNode::GetN(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputN(0);
}

const Sizes* TempVarNode::LocalM(unsigned int num) const
{
  if (num > 0)
    throw;
  return m_mlsize;
}

const Sizes* TempVarNode::LocalN(unsigned int num) const
{
  if (num > 0)
    throw;
  return m_nlsize;
}
#elif DOTENSORS
const Dim TempVarNode::NumDims(unsigned int num) const
{
  return InputNumDims(0);
}

const Sizes* TempVarNode::Len(unsigned int num, Dim dim) const
{
  if (num > 0)
    throw;
  return InputLen(0, dim);
}

const Sizes* TempVarNode::LocalLen(unsigned int num, Dim dim) const
{
  return &(m_lsizes[dim]);
}
#endif

Name TempVarNode::GetName(unsigned int num) const
{
  if (num > 0)
    throw; 
  Name tmp;
  if (m_name.empty()) {
    tmp = GetInputName(0);
    tmp.m_name = tmp.m_name + "temp";
  }
  else {
    tmp.m_name = m_name;
  }
#if DODM
  tmp.m_type = m_distType;
#endif
  return tmp;
}

void TempVarNode::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
#if DODM
  WRITE(m_distType);
#endif
  out << m_name << endl;
}

void TempVarNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLANode::UnflattenCore(in, info);
#if DODM
  READ(m_distType);
#endif
  getline(in,m_name);
}

void TempVarNode::ClearSizeCache()
{
#if TWOD
  if (!m_mlsize)
    return;
  delete m_mlsize;
  m_mlsize = NULL;
  delete m_nlsize;
  m_nlsize = NULL;
#elif DOTENSORS
  if (m_lsizes)
    return;
  delete [] m_lsizes;
  m_lsizes = NULL;
#endif
}

void TempVarNode::BuildSizeCache()
{
#if TWOD
  if (m_mlsize)
    return;
  m_mlsize = new Sizes;
  m_nlsize = new Sizes;
#if DODM
  GetLocalSizes(m_distType, GetM(0), GetN(0), *m_mlsize, *m_nlsize);
#elif DOBLIS
  *m_mlsize = *GetM(0);
  *m_nlsize = *GetN(0);
#else
sdlkfj
#endif
#elif DOTENSORS

  if (m_lsizes)
    return;
 Dim numDims = InputNumDims(0);
 m_lsizes = new Sizes[numDims];
 
 for (Dim dim = 0; dim < numDims; ++dim)
   GetLocalSizes(m_distType, dim, InputLen(0,dim), m_lsizes+dim);
 
#endif
}

#if TWOD
#if DOELEM
void MakeTrapNode::SanityCheck()
{
  DLAOp<1,1>::SanityCheck();
  if (m_inputs.size() != 1)
    cout << "m_inputs.size() != 1\n";
}

void MakeTrapNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<1,1>::Prop();
    m_cost = 0;
  }
}

void MakeTrapNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const MakeTrapNode *node = (MakeTrapNode*)orig;
  m_side = node->m_side;
  m_tri = node->m_tri;
  m_offset = node->m_offset;
  DLANode::Duplicate(orig,shallow,possMerging);
}

void MakeTrapNode::PrintCode(IndStream &out)
{  
  out.Indent();
  *out << "MakeTrapezoidal( " 
  << SideToStr(m_side) << ", "
  << TriToStr(m_tri) << ", " 
  << m_offset << ", " 
  << GetNameStr(0)
  << " );\n";
}

bool MakeTrapNode::CanTrans() const
{
  return ((DLANode*)Input(0))->CanTrans();
}

void MakeTrapNode::FlattenCore(ofstream &out) const
{
  DLAOp<1,1>::FlattenCore(out);
  WRITE(m_side);
  WRITE(m_tri);
  WRITE(m_offset);
}


void MakeTrapNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<1,1>::UnflattenCore(in,info);
  READ(m_side);
  READ(m_tri);
  READ(m_offset);
}

bool MoveMakeTrap::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != MakeTrapNode::GetClass()) {
    return false;
  }
  //This only works if we have one child.
  //The idea is that this node is on the abstract algorithm, but in code it
  // only gets applied to the final data that gets input into computation, 
  // not the middle redistributions or the original MC_MR data
  if (node->m_children.size() != 1) {
    cout << "MakeTrapNode has more than one child!\n";
    throw;
  }
  if (node->Child(0)->GetNodeClass() == RedistNode::GetClass())
    return true;
  return false;
}

void MoveMakeTrap::Apply(Poss *poss, Node *node) const
{
  Node *child = node->Child(0);
  if (child->GetNodeClass() != RedistNode::GetClass())
    throw;
  node->RedirectChildren(node->Input(0), node->InputConnNum(0));
  child->RedirectChildren(node, 0);
  node->ChangeInput2Way(node->Input(0), node->InputConnNum(0), child, 0);
  if (node->m_inputs.size() != 1)
    throw;
}

#endif

#if DOBLIS||DOELEM
bool RemoveScaleByOne::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != ScaleNode::GetClass())
    return false;
  const ScaleNode *scale = (ScaleNode*)node;
  if (scale->m_val == COEFONE) {
    return true;
  }
  else
    return false;
}

void RemoveScaleByOne::Apply(Poss *poss, Node *node) const
{
  node->RedirectChildren(node->Input(0), node->InputConnNum(0));
  node->m_poss->DeleteChildAndCleanUp(node);
}

void ScaleTrapNode::SanityCheck()
{
  DLAOp<1,1>::SanityCheck();
  if (m_inputs.size() != 1)
    cout << "m_inputs.size() != 1\n";
}

void ScaleTrapNode::PrintCode(IndStream &out)
{  
  Layer layer = GetLayer();
  out.Indent();
  if (layer == DMLAYER) {
    *out << "ScaleTrapezoid( ";
    out << m_val;
    *out << ", "
	 << SideToStr(m_side) << ", "
	 << TriToStr(m_tri) << ", 0, "
	 << GetNameStr(0) << " );\n";
  }
  else if (layer == S1LAYER 
	   || layer == S2LAYER
	   || layer == S3LAYER) 
    {
      //      bli_obj_set_struc( BLIS_TRIANGULAR, L_10_1 );                                  
      //      bli_obj_set_uplo( BLIS_LOWER, L_10_1 );  
      string name = GetInputNameStr(0);
      *out << "bli_obj_set_struc( BLIS_TRIANGULAR, " << name << " );\n";
      out.Indent();
      *out << "bli_obj_set_uplo( BLIS_LOWER, " << name << " );\n";
      out.Indent();
      *out << "bli_scalm( ";
      out << m_val;
      *out << ", &" << name << " );\n";
      out.Indent();
      *out << "bli_obj_set_struc( BLIS_GENERAL, " << name << " );\n";
      out.Indent();
      *out << "bli_obj_set_uplo( BLIS_DENSE, " << name << " );\n";
    }
}

void ScaleTrapNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<1,1>::Duplicate(orig, shallow, possMerging);
  ScaleTrapNode *scal = (ScaleTrapNode*)orig;
  m_val = scal->m_val;
  m_side = scal->m_side;
  m_tri = scal->m_tri;
}

#if DOELEM
const DistType& ScaleTrapNode::GetDistType(unsigned int num) const
{
  if (num > 0)
    throw;
  return ((DLANode*)Input(0))->GetDistType(InputConnNum(0));
}
#endif


void ScaleTrapNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<1,1>::Prop();
    m_cost = 0;
  }
}

void ScaleTrapNode::FlattenCore(ofstream &out) const
{
  DLAOp<1,1>::FlattenCore(out);
  WRITE(m_side);
  WRITE(m_tri);
  WRITE(m_val);
}

void ScaleTrapNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<1,1>::UnflattenCore(in, info);
  READ(m_side);
  READ(m_tri);
  READ(m_val);
}

ScaleNode::ScaleNode(Layer layer, Coef val) 
  : m_val(val) 
{
  SetLayer(layer); 
}

void ScaleNode::SanityCheck()
{
  DLAOp<1,1>::SanityCheck();
  if (m_inputs.size() != 1)
    cout << "m_inputs.size() != 1\n";
}

void ScaleNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<1,1>::Prop();
    m_cost = 0;
  }
}

void ScaleNode::PrintCode(IndStream &out)
{  
  out.Indent();
  if (GetLayer() == DMLAYER || GetLayer() == ABSLAYER) {
    *out << "Scale( ";
  }
  else if (GetLayer() == S1LAYER || GetLayer() == S2LAYER || GetLayer() == S3LAYER)
    *out << "bli_scalm( ";

  out << m_val;
  *out << ", "
	 << GetNameStr(0) << " );\n";

}

void ScaleNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<1,1>::Duplicate(orig, shallow, possMerging);
  ScaleNode *scal = (ScaleNode*)orig;
  m_val = scal->m_val;
}

#if DODM
const DistType& ScaleNode::GetDistType(unsigned int num) const
{
  if (num > 0)
    throw;
  return ((DLANode*)Input(0))->GetDistType(InputConnNum(0));
}
#endif

void ScaleNode::FlattenCore(ofstream &out) const
{
  DLAOp<1,1>::FlattenCore(out);
  WRITE(m_val);
}


void ScaleNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<1,1>::UnflattenCore(in, info);
  READ(m_val);
}
#endif


#if DOELEM
ViewPan::~ViewPan()
{
  if (m_sizes) {
    delete m_sizes;
    m_sizes = NULL;
    delete m_lsizes;
    m_lsizes = NULL;
  }
}

void ViewPan::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (m_inputs.size() != 2)
      throw;
    Input(0)->Prop();
    Input(1)->Prop();
    
    if (m_isVert) {
      if (*GetInputN(0) != *GetInputN(1))
        throw;
    }
    else {
      if (*GetInputM(0) != *GetInputM(1)) {
        throw;
      }
    }
  }
}

void ViewPan::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  m_isVert = ((ViewPan*)orig)->m_isVert;
  m_name = ((ViewPan*)orig)->m_name;
}

void ViewPan::ClearSizeCache()
{
  if (m_sizes) {
    delete m_sizes;
    m_sizes = NULL;
    delete m_lsizes;
    m_lsizes = NULL;
  }
}

void ViewPan::BuildSizeCache()
{
  if (m_sizes)
    return;
  m_sizes = new Sizes;
  m_lsizes = new Sizes;
  DistType t = GetDistType(0);
  if (m_isVert) {
    const Sizes *size1 = GetInputM(0);
    const Sizes *size2 = GetInputM(1);
    const Sizes *nSizes = GetInputN(0);
    Sizes tmp;
    m_sizes->PairwiseSum(*size1,*size2);
    GetLocalSizes(t, m_sizes, nSizes, *m_lsizes, tmp);
  }
  else {
    const Sizes *size1 = GetInputN(0);
    const Sizes *size2 = GetInputN(1);
    if (size1->NumSizes() != size2->NumSizes())
      throw;
    const Sizes *mSizes = GetInputM(0);
    m_sizes->PairwiseSum(*size1,*size2);
    Sizes tmp;
    GetLocalSizes(t, mSizes, m_sizes, tmp, *m_lsizes);
  }
}

const Sizes* ViewPan::GetM(unsigned int num) const
{
  if (num >  0) 
    throw;
  if (m_isVert)
    return m_sizes;
  else
    return GetInputM(0);
}

const Sizes* ViewPan::GetN(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_isVert)
    return GetInputN(0);
  else
    return m_sizes;
}


const Sizes* ViewPan::LocalM(unsigned int num) const
{
  if (num >  0) 
    throw;
  if (m_isVert)
    return m_lsizes;
  else
    return InputLocalM(0);
}

const Sizes* ViewPan::LocalN(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_isVert)
    return InputLocalN(0);
  else
    return m_lsizes;
}

Name ViewPan::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  Name name = GetInputName(0);
  name.m_name = m_name;
  return name;
}

void ViewPan::PrintCode(IndStream &out)
{
  out.Indent();
  *out << GetNameStr(0);
  if (m_isVert)
    *out << ".View2x1 ( "
    << GetInputName(0).str() << ",\n"
    << out.Tabs(1)
    << GetInputName(1).str() << " );\n";
  else
    *out << ".View1x2 ( "
    << GetInputName(0).str() << ", "
    << GetInputName(1).str() << " );\n";
}

void ViewPan::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  WRITE(m_isVert);
  out << m_name << endl;
}


void ViewPan::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLANode::UnflattenCore(in, info);
  READ(m_isVert);
  getline(in, m_name);
}

ViewAroundDiag::~ViewAroundDiag()
{
  if (m_sizes0) {
    delete m_sizes0;
    m_sizes0 = NULL;
    delete m_lsizes0;
    m_lsizes0 = NULL;
    delete m_sizes1;
    m_sizes1 = NULL;
    delete m_lsizes1;
    m_lsizes1 = NULL;
  }
}

void ViewAroundDiag::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (m_inputs.size() != 3)
      throw;
    Input(0)->Prop();
    Input(1)->Prop();
    Input(2)->Prop();
    
    if (m_isVert) {
      if (*GetInputN(0) != *GetInputN(1))
        throw;
      if (*GetInputN(1) != *GetInputN(2))
        throw;
    }
    else {
      if (*GetInputM(0) != *GetInputM(1))
        throw;
      if (*GetInputM(1) != *GetInputM(2))
        throw;
    }
  }
}

void ViewAroundDiag::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  m_isVert = ((ViewAroundDiag*)orig)->m_isVert;
  m_name = ((ViewAroundDiag*)orig)->m_name;
}

void ViewAroundDiag::ClearSizeCache()
{
  if (m_sizes0) {
    delete m_sizes0;
    m_sizes0 = NULL;
    delete m_lsizes0;
    m_lsizes0 = NULL;
    delete m_sizes1;
    m_sizes1 = NULL;
    delete m_lsizes1;
    m_lsizes1 = NULL;
  }
}

void ViewAroundDiag::BuildSizeCache()
{
  if (m_sizes0)
    return;
  m_sizes0 = new Sizes;
  m_lsizes0 = new Sizes;
  m_sizes1 = new Sizes;
  m_lsizes1 = new Sizes;
  DistType t = GetDistType(0);
  if (m_isVert) {
    const Sizes *size0 = GetInputM(0);
    const Sizes *size1 = GetInputM(1);
    const Sizes *size2 = GetInputM(2);
    const Sizes *nSizes = GetInputN(0);
    m_sizes0->PairwiseSum(*size0, *size1);
    m_sizes1->PairwiseSum(*size1, *size2);
    Sizes tmp;
    ::GetLocalSizes(t, m_sizes0, nSizes, *m_lsizes0, tmp);
    ::GetLocalSizes(t, m_sizes1, nSizes, *m_lsizes1, tmp);
  }
  else {
    const Sizes *size0 = GetInputN(0);
    const Sizes *size1 = GetInputN(1);
    const Sizes *size2 = GetInputN(2);
    const Sizes *mSizes = GetInputM(0);
    m_sizes0->PairwiseSum(*size0, *size1);
    m_sizes1->PairwiseSum(*size1, *size2);
    Sizes tmp;
    ::GetLocalSizes(t, mSizes, m_sizes0, tmp, *m_lsizes0);
    ::GetLocalSizes(t, mSizes, m_sizes1, tmp, *m_lsizes1);
  }
}

const Sizes* ViewAroundDiag::GetM(unsigned int num) const
{
  if (num >  1) 
    throw;
  if (m_isVert) {
    if (num == 0)
      return m_sizes0;
    else
      return m_sizes1;
  }
  else
    return GetInputM(0);
}

const Sizes* ViewAroundDiag::GetN(unsigned int num) const
{
  if (num > 1)
    throw;
  if (m_isVert)
    return GetInputN(0);
  else {
    if (num == 0)
      return m_sizes0;
    else
      return m_sizes1;
  }
}


const Sizes* ViewAroundDiag::LocalM(unsigned int num) const
{
  if (num >  1) 
    throw;
  if (m_isVert) {
    if (num == 0)
      return m_lsizes0;
    else
      return m_lsizes1;
  }
  else
    return InputLocalM(0);
}

const Sizes* ViewAroundDiag::LocalN(unsigned int num) const
{
  if (num > 1)
    throw;
  if (m_isVert)
    return InputLocalN(0);
  else {
    if (num == 0)
      return m_lsizes0;
    else
      return m_lsizes1;
  }
}

Name ViewAroundDiag::GetName(unsigned int num) const
{
  if (num > 1)
    throw;
  Name name = GetInputName(0);
  name.m_name = m_name;
  if (m_isVert) {
    if (num == 0) {
      name.m_name += "Above";
    }
    else {
      name.m_name += "Below";
    }
  }
  else {
    if (num == 0) {
      name.m_name += "Left";
    }
    else {
      name.m_name += "Right";
    }
  }
  return name;
}

void ViewAroundDiag::PrintCode(IndStream &out)
{
  if (m_isVert) {
    out.Indent();
    *out << GetNameStr(0)
    << ".View2x1 ( "
    << GetInputName(0).str() << ",\n"
    << out.Tabs(1)
    << GetInputName(1).str() << " );\n";
    out.Indent();
    *out << GetNameStr(1)
    << ".View2x1 ( "
    << GetInputName(1).str() << ",\n"
    << out.Tabs(1)
    << GetInputName(2).str() << " );\n";
  }
  else {
    out.Indent();
    *out << GetNameStr(0)
    << ".View1x2 ( "
    << GetInputName(0).str() << ", "
    << GetInputName(1).str() << " );\n";
    out.Indent();
    *out << GetNameStr(1)
    << ".View1x2 ( "
    << GetInputName(1).str() << ", "
    << GetInputName(2).str() << " );\n";
  }
}

void ViewAroundDiag::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  WRITE(m_isVert);
  out << m_name << endl;
}

void ViewAroundDiag::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLANode::UnflattenCore(in, info);
  READ(m_isVert);
  getline(in, m_name);
}

void ViewAroundDiagCombine::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<5,3>::Prop();
    m_cost = 0;
  }
}

void ViewAroundDiagCombine::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<5,3>::Duplicate(orig, shallow, possMerging);
  m_isVert = ((ViewAroundDiagCombine*)orig)->m_isVert;
}

void ViewAroundDiagCombine::FlattenCore(ofstream &out) const
{
  DLAOp<5,3>::FlattenCore(out);
  WRITE(m_isVert);
}


void ViewAroundDiagCombine::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<5,3>::UnflattenCore(in, info);
  READ(m_isVert);
}
#endif

Name ViewTL::GetName(unsigned int num) const
{
  Name name = GetInputName(0);
  name.m_name += "TL";
  return name;    
}

 void ViewTL::Prop()
{
  if (!IsValidCost(m_cost)) {
    Input(0)->Prop();
    Input(1)->Prop();
    Input(2)->Prop();
    //    DLANode::Prop();
    m_cost = ZERO;
  }
}

void ViewTL::SanityCheck()
{
  DLANode::SanityCheck();
}

void ViewTL::PrintCode(IndStream &out)
{
  if (GetLayer() == S3LAYER) {
    string name = GetNameStr(0);
    out.Indent();
    *out << "obj_t " << name << ", " << name << "tmp;\n";
    out.Indent();
    *out << "bli_acquire_mpart_l2r( BLIS_SUBPART1, 0, bli_obj_width( " 
	 << GetInputNameStr(1) << " ), &" << GetInputNameStr(0)
	 << ", &" << name << "tmp );\n";
    out.Indent();
    *out << "bli_acquire_mpart_t2b( BLIS_SUBPART1, 0, bli_obj_length( " 
	 << GetInputNameStr(2) << " ), &" << name
	 << "tmp, &" << name << " );\n";
  }
  else
    throw;
}


void ViewTLCombine::Prop()
{
  if (!IsValidCost(m_cost)) {
    m_cost = ZERO;
    Input(0)->Prop();
    Input(1)->Prop();
  }
}
#endif
