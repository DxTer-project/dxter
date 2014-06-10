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


#include "layers.h"
#if DOELEM||DOBLIS
#include "base.h"
#include "blas.h"
#include "elemRedist.h"
#include "string.h"
#include "helperNodes.h"

using namespace std;

Axpy::Axpy(Layer layer, Coef coeff) 
: m_coeff(coeff)
{ 
  SetLayer(layer); 
#if DOBLIS
  m_comm = CORECOMM;
#endif
}

NodeType Axpy::GetType() const
{
  return "Axpy " + LayerNumToStr(GetLayer()) 
#if DOBLIS
    + CommToStr(m_comm);
#else
  ;
#endif
}

void Axpy::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<2,1>::Duplicate(orig, shallow, possMerging);
  const Axpy *axpy = (Axpy*)orig;
  m_coeff = axpy->m_coeff;
#if DOBLIS
  m_comm = axpy->m_comm;
#endif
}

void Axpy::FlattenCore(ofstream &out) const
{
  DLAOp<2,1>::FlattenCore(out);
  WRITE(m_coeff);
#if DOBLIS
  WRITE(m_comm);
#endif
}

void Axpy::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<2,1>::UnflattenCore(in, info);
  READ(m_coeff);
#if DOBLIS
  READ(m_comm);
#endif
}

#if DOELEM
const DistType& Axpy::GetDistType(unsigned int num) const 
{ 
  if (m_layer == DMLAYER || m_layer == ABSLAYER)
    return MC_MR; 
  else if (m_layer == SMLAYER)
    return InputDistType(1);
  else
    throw;
}
#endif

Phase Axpy::MaxPhase() const 
{
#if DODPPHASE
  if (m_layer == DMLAYER)
    return DPPHASE;
  else if (m_layer == SMLAYER)
    return NUMPHASES;
  else
    throw;
#elif DOSM||DOSQM
  if (m_layer == ABSLAYER)
    return SR1PHASE;
  else if (m_layer == S3LAYER)
    return NUMPHASES;
  else
    throw;
#else
  throw;
#endif
}

bool Axpy::ShouldCullDP() const 
{
#if DODPPHASE
  return m_layer == DMLAYER;
#else
  throw;
#endif
}

bool Axpy::DoNotCullDP() const 
{
#if DODPPHASE
  return m_layer == DMLAYER;
#else
  throw;
#endif
}


void Axpy::SanityCheck()
{
  DLAOp<2,1>::SanityCheck();
  if (m_inputs.size() != 2)
    cout << "2 m_inputs.size() != 2\n";
#if DODPPHASE
  else if (m_layer == DMLAYER) {
    if (InputDistType(0) != D_MC_MR)
      throw;
    if (InputDistType(1) != D_MC_MR)
      throw;
  }
  else if (m_layer == SMLAYER) {
    if (InputDistType(0) != InputDistType(1))
      m_poss->MarkInsane();
  }
#if DOBLIS
  if (m_comm != CORECOMM)
    throw;
#endif
#else

#endif
}

void Axpy::PrintCode(IndStream &out)
{
  out.Indent();
  switch (GetLayer()) {
  case (SMLAYER):
    *out << "Axpy( ";
    out << m_coeff;
    *out << ", "
	 << GetInputName(0).str() << ","
	 << "\n" << out.Tabs(2) << GetInputName(1).str()
	 << ");\n";
    break;
#if DOBLIS
  case (S3LAYER):
    if (m_comm == CORECOMM) {
      *out << "bli_axpym( ";
    }
    else {
      *out << "reduce( " << CommToStr(m_comm) << ", ";
    }
    out << m_coeff;
    *out << ", &"
	 << GetInputName(0).str() << ","
	 << "\n" << out.Tabs(2) << "&" << GetInputName(1).str()
	 << ");\n";
    break;
#endif
  default:
    throw;
  }
}

void Axpy::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();
    if (m_layer == DMLAYER)
      m_cost = 0;
    else if (m_layer == SMLAYER)
      m_cost = GetCost(SMLAYER, LocalM(0), LocalN(0));
#if DOBLIS
    else if (m_layer == S3LAYER) {
      m_cost = GetCost(S3LAYER, LocalM(0), LocalN(0));
      if (m_comm != CORECOMM) {
	m_cost /= NumCoresInComm(m_comm);
      }
    }
#endif
    else if (m_layer == ABSLAYER)
      m_cost = ZERO;
    else {
      cout << LayerNumToStr(m_layer) << endl;
      throw;
    }
  }
}

Cost Axpy::GetCost(Layer layer, const Sizes *localMs, const Sizes *localNs)
{
  if (layer == SMLAYER)
    return TWO * GAMMA * localMs->SumProds11(*localNs);
  else if (layer == S3LAYER)
    return TWO * GAMMA * localMs->SumProds11(*localNs) + (2*PSIWVAL+PSIRVAL) * localNs->Sum();
  else
    throw;
}

#if DOELEM
bool FindDistUp(const RedistNode *node, DistType type)
{
  if (node->m_destType == type)
    return true;
  Node *input = node->Input(0);
  if (input->GetNodeClass() == RedistNode::GetClass())
    return FindDistUp((RedistNode*)input, type);
  else
    return node->InputDistType(0) == type;
}

bool FindDistDown(const Node *node, DistType type)
{
  if (node->GetNodeClass() != RedistNode::GetClass()) 
    return false;
  if (((RedistNode*)node)->m_destType == type)
    return true;
  else {
    NodeConnVecConstIter iter = node->m_children.begin();
    for(; iter != node->m_children.end(); ++iter) {
      const Node *child = (*iter)->m_n;
      if (FindDistDown((RedistNode*)child,type))
	return true;
    }
  }
  return false;
}


bool CheckInput(const Node *node, unsigned int inNum, DistType type, bool skipFirstRedist) 
{
  const Node *in = node->Input(inNum);
  if (skipFirstRedist && in->GetNodeClass() == RedistNode::GetClass()) {
    return CheckInput(in, 0, type, false);
  }
  DistType inType = ((DLANode*)node)->InputDistType(inNum);  
  if (inType == type)
    return true;
  if (in->GetNodeClass() == RedistNode::GetClass()) {
    return FindDistUp((RedistNode*)in, type);
  }
  return false;
}
#endif

#if DOELEM
bool DistAxpyToLocalAxpy::WorthApplying(const Node *node) const
{
  if (node->GetNodeClass() != Axpy::GetClass())
    return false;
  if (((Axpy*)node)->GetLayer() != DMLAYER)
    return false;
  if (CheckInput(node, 0, m_type, true))
    return true;
  if (CheckInput(node, 1, m_type, true))
    return true;
  NodeConnVecConstIter iter = node->m_children.begin();
  for(; iter != node->m_children.end(); ++iter) {
    const Node *child = (*iter)->m_n;
    if (FindDistDown((RedistNode*)child,m_type))
      return true;
  }
  return false;
}

bool DistAxpyToLocalAxpy::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != Axpy::GetClass())
    return false;
  if (((Axpy*)node)->GetLayer() != DMLAYER)
    return false;
  return true;
}

void DistAxpyToLocalAxpy::Apply(Poss *poss, Node *node) const
{
  RedistNode *node1 = new RedistNode(m_type);
  RedistNode *node2 = new RedistNode(m_type);
  Axpy *node3 = new Axpy(SMLAYER, (((Axpy*)node)->m_coeff));
  RedistNode *node4 = new RedistNode(D_MC_MR);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node->Input(1),node->InputConnNum(1));
  node3->AddInput(node1,0);
  node3->AddInput(node2,0);
  node4->AddInput(node3,0);
  poss->AddNodes(4, node1, node2, node3, node4);
  node->RedirectChildren(node4,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

Cost DistAxpyToLocalAxpy::RHSCostEstimate(const Node *node) const
{
  const Axpy *axpy = (Axpy*)node;
  const Sizes *A1 = axpy->GetInputM(0);
  const Sizes *A2 = axpy->GetInputN(0);
  Cost cost = 0;
  if (m_type != D_MC_MR) {
    cost += 2 * RedistNode::GetCost(D_MC_MR, m_type, A1, A2);
    cost += RedistNode::GetCost(m_type, D_MC_MR, A1, A2);
  }
  Sizes localM;
  Sizes localN;
  GetLocalSizes(m_type, A1, A2, localM, localN);
  cost += Axpy::GetCost(SMLAYER, &localM, &localN);

  return cost;
}
#endif

bool AxpyToBLASAxpy::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != Axpy::GetClass())
    return false;
  return ((DLANode*)node)->GetLayer() == m_fromLayer;
}

void AxpyToBLASAxpy::Apply(Poss *poss, Node *node) const
{
  ((DLANode*)node)->SetLayer(m_toLayer);
}

string AxpyLowerLayer::GetType() const 
{ 
  return "Axpy lower layer " + LayerNumToStr(m_fromLayer) 
    + " to " + LayerNumToStr(m_toLayer);
}

bool AxpyLowerLayer::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Axpy::GetClass()) {
    const Axpy *axpy = (Axpy*)node;
    return (axpy->GetLayer() == m_fromLayer);
  }
  return false;
  
}

void AxpyLowerLayer::Apply(Poss *poss, Node *node) const
{
  Axpy *axpy = (Axpy*)node;
  axpy->SetLayer(m_toLayer);
}


NodeType Scal::GetType() const 
{
  return "Scal" + LayerNumToStr(GetLayer());
}

#if DOELEM
const DistType& Scal::GetDistType(unsigned int num) const
{ 
  if (m_layer == DMLAYER)
    return MC_MR; 
  else if (m_layer == SMLAYER)
    return InputDistType(1);
  else if (m_layer == S3LAYER)
    return InputDistType(1);
  else
    throw;
}
#endif

Phase Scal::MaxPhase() const 
{
#if DODPPHASE
  if (m_layer == DMLAYER)
    return DPPHASE;
  else if (m_layer == SMLAYER)
    return NUMPHASES;
  else
    throw;
#else
  throw;
#endif
}

bool Scal::ShouldCullDP() const 
{
#if DODPPHASE
  return m_layer == DMLAYER;
#else
  throw;
#endif
}

void Scal::SanityCheck()
{
  DLAOp<2,1>::SanityCheck();
  if (m_inputs.size() != 2) {
    cout << "1 m_inputs.size() != 2\n";
    throw;
  }
  if (!((DLANode*)Input(0))->IsScalar(InputConnNum(0)))
    throw;
#if DODPPHASE
  else if (m_layer == DMLAYER) {
    if (InputDistType(0) != D_STAR_STAR)
      throw;
    if (InputDistType(1) != D_MC_MR)
      throw;
  }
  else if (m_layer == SMLAYER) {
    if (InputDistType(0) != D_STAR_STAR)
      cout << "input not D_STAR_STAR";
  }
  else
    throw;
#else
  throw;
#endif
}

void Scal::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "Scal( " << GetInputName(0).str() << ".LocalBuffer(0,0), " << GetInputName(1).str() << ");\n";
}

void Scal::Prop()
{
  if (!IsValidCost(m_cost)) {
    Scal::Prop();
    if (m_layer == SMLAYER)
      m_cost = GAMMA * LocalM(1)->SumProds11(*LocalN(1));
    else if (m_layer == DMLAYER)
      m_cost = 0;
    else
      throw;
  }
}

#if DOELEM
string DistScalToLocalScal::GetType() const
{
  return "Distributed Scal to Local Scal " + DistTypeToStr(m_type);
}

bool DistScalToLocalScal::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != Scal::GetClass())
    return false;
  return (((Scal*)node)->GetLayer() == DMLAYER);
}

void DistScalToLocalScal::Apply(Poss *poss, Node *node) const
{
  RedistNode *node1 = new RedistNode(D_STAR_STAR);
  RedistNode *node2 = new RedistNode(m_type);
  Scal *node3 = new Scal(SMLAYER);
  RedistNode *node4 = new RedistNode(D_MC_MR);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node->Input(1),node->InputConnNum(1));
  node3->AddInput(node1,0);
  node3->AddInput(node2,0);
  node4->AddInput(node3,0);
  poss->AddNodes(4, node1, node2, node3, node4);
  node->RedirectChildren(node4,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif


Phase ConstScal::MaxPhase() const 
{
#if DODPPHASE
  if (m_layer == DMLAYER)
    return DPPHASE;
  else if (m_layer == SMLAYER)
    return NUMPHASES;
  else
    throw;
#else
  throw;
#endif
}

#if DOELEM
const DistType& ConstScal::GetDistType(unsigned int num) const
{ 
  if (m_layer == DMLAYER)
    return MC_MR; 
  else if (m_layer == SMLAYER)
    return InputDistType(0);
  else
    throw;
}
#endif

bool ConstScal::ShouldCullDP() const 
{
#if DODPPHASE
  return m_layer == DMLAYER;
#else
  throw;
#endif
}


void ConstScal::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<1,1>::Duplicate(orig, shallow, possMerging);
  const ConstScal *scal = (ConstScal*)orig;
  m_alpha = scal->m_alpha;
}

NodeType ConstScal::GetType() const 
{
  return "ConstScal" + LayerNumToStr(GetLayer());
}

void ConstScal::FlattenCore(ofstream &out) const
{
  DLAOp<1,1>::FlattenCore(out);
  WRITE(m_alpha);
}


void ConstScal::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<1,1>::UnflattenCore(in, info);
  READ(m_alpha);
}

void ConstScal::SanityCheck()
{
  DLAOp<1,1>::SanityCheck();
  if (m_inputs.size() != 1) {
    cout << "1 m_inputs.size() != 1\n";
    throw;
  }
#if DOELEM
  else if (m_layer == DMLAYER) {
    if (InputDistType(1) != D_MC_MR)
      throw;
  }
#endif
}

void ConstScal::PrintCode(IndStream &out)
{
  out.Indent();
  if (m_layer == DMLAYER) {
    *out << "DistConstScal( ";
    out << m_alpha;
    *out << ", " << GetInputName(1).str() << ");\n";
  }
  else if (m_layer == SMLAYER) {
    *out << "ConstScal( ";
    out << m_alpha;
    *out << ", " << GetInputName(0).str() << ");\n";
  }
  else
    throw;
}

void ConstScal::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<1,1>::Prop();
    if (m_layer == SMLAYER)
      m_cost = GAMMA * LocalM(0)->SumProds11(*LocalN(0));
    else if (m_layer == DMLAYER)
      m_cost = 0;
    else
      throw;
  }
}

#if DODPHASE
string DistConstScalToLocalConstScal::GetType() const
{
  return "Distributed ConstScal to Local ConstScal " + DistTypeToStr(m_type);
}

bool DistConstScalToLocalConstScal::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != DistConstScal::GetClass())
    return false;
  return true;
}

void DistConstScalToLocalConstScal::Apply(Poss *poss, Node *node) const
{
  DistConstScal *scal = (DistConstScal*)node;
  RedistNode *node2 = new RedistNode(m_type);
  LocalConstScal *node3 = new LocalConstScal(scal->m_alpha);
  RedistNode *node4 = new RedistNode(D_MC_MR);
  node2->AddInput(node->Input(1),node->InputConnNum(1));
  node3->AddInput(node2,0);
  node4->AddInput(node3,0);
  poss->AddNodes(3, node2, node3, node4);
  node->RedirectChildren(node4,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif
#endif
