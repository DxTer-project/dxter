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
#include "contraction.h"

#if DOTENSORS
#include "tensorRedist.h"
#include "loopSupport.h"
#include "helperNodes.h"
#include "yxpby.h"

typedef vector<unsigned int *> ContType;
typedef ContType::iterator ContTypeIter;

void MatchDistsAndFillInWithStar(string indices, 
				 const DistType &matchingDists, string matchingIndices,
				 DistType &final);
/*
void RecursivelyFindDistributions(DimVec *dists, Dim thisDim, 
				  const DistType &AType, const DimVec &ADims,
				  const DistType &BType, const DimVec &BDims,
				  DimSet &usedDims, Dim numDims,
				  ContType *distOptions);
void AddUnusedDimsForDistType(DimVec *dists, unsigned int *distEntries,
			      Dim numIndices, 
			       DimSet &unUsedDims,
			       ContType *distOptions);

void MatchDistsAndFillIn(string indices, 
			 const DistType &matchingDists, string matchingIndices,
			 unsigned int *fillInDists, const DimVec &fillInDims,
			 DistType &final);
*/

Contraction::Contraction(Layer layer, Coef alpha, Coef beta, Type type, 
	      string AIndices, string BIndices, string CIndices,
	      string contIndices)
:
  m_alpha(alpha),
  m_beta (beta),
  m_type (type),
  m_AIndices(AIndices),
  m_BIndices(BIndices),
  m_CIndices(CIndices),
  m_contIndices(contIndices)
{
  SetLayer(layer);
}

Node* Contraction::BlankInst()
{
  return new Contraction(ABSLAYER, COEFONE, COEFONE, REAL, "", "", "", "");
}

NodeType Contraction::GetType() const
{
  return "Contraction " + m_AIndices + " " + m_BIndices + " " 
    + m_CIndices + " " + m_contIndices
    + LayerNumToStr(GetLayer());
}

void Contraction::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const Contraction *cont = (Contraction*)orig;
  m_alpha = cont->m_alpha;
  m_beta = cont->m_beta;
  m_type = cont->m_type;
  m_AIndices = cont->m_AIndices;
  m_BIndices = cont->m_BIndices;
  m_CIndices = cont->m_CIndices;
  m_contIndices = cont->m_contIndices;
}

void Contraction::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_alpha);
  WRITE(m_beta);
  WRITE(m_type);
  out << m_AIndices <<endl;
  out << m_BIndices <<endl;
  out << m_CIndices <<endl;
  out << m_contIndices << endl;
}

void Contraction::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  READ(m_alpha);
  READ(m_beta);
  READ(m_type);
  getline(in, m_AIndices);
  getline(in, m_BIndices);
  getline(in, m_CIndices);
  getline(in, m_contIndices);
}

void Contraction::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();
    m_cost = 0;

    if(InputNumDims(0) != m_AIndices.size())
      throw;

    if(InputNumDims(1) != m_BIndices.size())
      throw;

    if(InputNumDims(2) != m_CIndices.size()) {
      cout << "num " << InputNumDims(2) << endl;
      cout << "dist " << InputDataType(2).m_dist.PrettyStr() << endl;
      cout << m_CIndices << endl;
      cout << "Input " << Input(2)->GetNodeClass() << endl;
      throw;
    }

    //    const DistType &CType = InputDataType(2).m_dist;

    string::iterator strIter = m_contIndices.begin();
    for(; strIter != m_contIndices.end(); ++strIter) {
      const char index = *strIter;
      if (m_AIndices.find(index) == string::npos)
	throw;
      if (m_BIndices.find(index) == string::npos)
	throw;
      //      if (m_CIndices.find(index) != string::npos)
      //	throw;
    }

    CheckInputTypesAlign();
    
    //    cout << "improve Contraction::Prop code\n";
    //    cout << "reflect in DistContToLocalContStatC::RHSCostEstimate\n";
    DimVec dims = MapIndicesToDims(m_contIndices,m_AIndices);
    const Sizes *sizes = InputLocalLen(2,0);
    Dim numDims = InputNumDims(2);
    unsigned int totNumIters = sizes->NumSizes();

    if (m_layer == ABSLAYER || m_layer == DMLAYER) {
      for(unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
	Cost temp = 2;
	for (Dim dim = 0; dim < numDims; ++dim) {
	  temp *= (*InputLen(2,dim))[iteration];
	}
	DimVecConstIter iter = dims.begin();
	for(; iter != dims.end(); ++iter) {
	  temp *= (*InputLen(0,*iter))[iteration];
	}
	m_cost += temp;
      }
    }
    else if (m_layer == SMLAYER) {
      for(unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
	Cost temp = 2;
	for (Dim dim = 0; dim < numDims; ++dim) {
	  temp *= (*InputLocalLen(2,dim))[iteration];
	}
	DimVecConstIter iter = dims.begin();
	for(; iter != dims.end(); ++iter) {
	  temp *= (*InputLocalLen(0,*iter))[iteration];
	}
	m_cost += temp;
      }
    }
    else
      throw;
  }
}

void Contraction::AddVariables(VarSet &set) const
{
  //Reflect in PrintCode
  DLANode::AddVariables(set);
  
  Var varA(IndexArrayType, m_AIndices);
  set.insert(varA);
  Var varB(IndexArrayType, m_BIndices);
  set.insert(varB);
  Var varC(IndexArrayType, m_CIndices);
  set.insert(varC);
}

void Contraction::CheckInputTypesAlign() const
{
  if (m_layer != SMLAYER)
    return;

  const DistType &AType = InputDataType(0).m_dist;
  const DistType &BType = InputDataType(1).m_dist;
  const DistType &CType = InputDataType(2).m_dist;
    
  Dim dim = 0;
  string::const_iterator strIter = m_CIndices.begin();
  for(; strIter != m_CIndices.end(); ++strIter, ++dim) {
    const char index = *strIter;
    size_t aFind = m_AIndices.find(index);
    if (aFind == string::npos) {
      size_t bFind = m_BIndices.find(index);
      if (bFind == string::npos) {
	if (!InputLen(2,dim)->AllOnes())
	  throw;
	//else we're contracting to a scalar or vector or matrix etc.
      }
      else {
	if (CType.m_dists[dim] != BType.m_dists[bFind]) {
	  cout << AType.PrettyStr() << endl;
	  cout << m_AIndices << endl;
	  cout << BType.PrettyStr() << endl;
	  cout << m_BIndices << endl;
	  cout << CType.PrettyStr() << endl;
	  cout << m_CIndices << endl;
	  cout << m_contIndices << endl;


	  cout << *strIter << " " << m_BIndices[bFind] << endl;
	  cout << CType.m_dists[dim].str() << endl;
	  cout << BType.m_dists[bFind].str() << endl;
	  throw;
	}
	if (m_contIndices.find(index) == string::npos &&
	    *InputLocalLen(2,dim) != *InputLocalLen(1,bFind)) {
	  cout << "C and B sizes don't align\n";
	  throw;
	}	  
      }
    }
    else {
      if (CType.m_dists[dim] != AType.m_dists[aFind]) {
	cout << "C\n";
        cout << CType.PrettyStr() << endl;
        cout << m_CIndices << endl;
        cout << dim << endl;
	cout << "A\n";
        cout << AType.PrettyStr() << endl;
        cout << aFind << endl;
        cout << m_AIndices << endl;
        cout << Input(0)->GetType() << endl;
        cout << Input(2)->GetType() << endl;
        throw;
      }
      if (m_contIndices.find(index) == string::npos && 
	  *InputLocalLen(2,dim) != *InputLocalLen(0,aFind)) {
	cout << "C and A sizes don't align\n";
	cout << "C:\n";
	InputLocalLen(2,dim)->Print();
	cout << "A:\n";
	InputLocalLen(0,aFind)->Print();
	cout << m_CIndices << endl;
	cout << m_AIndices << endl;
	cout << m_contIndices << endl;
	throw;
      }
      size_t bFind = m_BIndices.find(index);
      if (bFind != string::npos) {
	if (AType.m_dists[aFind].IsStar())
	  throw;
	if (BType.m_dists[bFind] != AType.m_dists[aFind])
	  throw;
	if (!InputLocalLen(2,dim)->AllOnes())
	  throw;
      }
    }
  }

  dim = 0;
  strIter = m_contIndices.begin();
  for(; strIter != m_contIndices.end(); ++strIter, ++dim) {
    const char index = *strIter;
    size_t aFind = m_AIndices.find(index);
    size_t bFind = m_BIndices.find(index);
    if (*InputLocalLen(0,aFind) != *InputLocalLen(1,bFind)) {
      cout << "A and B contraction sizes don't align\n";
    }
    if (AType.m_dists[aFind] != BType.m_dists[bFind]) {
      cout << "A and B contraction distributions don't align\n";
    }
  }
}


void Contraction::PrintCode(IndStream &out)
{
  Name in0 = GetInputName(0);
  Name in1 = GetInputName(1);
  Name in2 = GetInputName(2);

  out.Indent();
  *out << "   // ";
  out << m_alpha;
  *out << " * " << in0.PrettyStr() + "_" + m_AIndices
       << " * " << in1.PrettyStr() + "_" + m_BIndices
       << " + ";
  out << m_beta;
  *out << " * " << in2.PrettyStr() + "_" + m_CIndices
       << endl;

  out.Indent();

  if (m_CIndices.find(m_contIndices[0]) == string::npos) 
    *out << "LocalContractAndLocalEliminate(";
  else
    *out << "LocalContract(";
  out << m_alpha;
  *out << ", " << in0.str()
       << ".LockedTensor(), " 
       << IndexArrayVarName(m_AIndices) << ",\n";
  out.Indent(1);
  *out << in1.str()
       << ".LockedTensor(), "
       << IndexArrayVarName(m_BIndices) << ",\n";
  out.Indent(1);
  out << m_beta;
  *out << ", " << in2.str() << ".Tensor(), "
       << IndexArrayVarName(m_CIndices) << ");\n";

  
}


Phase Contraction::MaxPhase() const
{
  switch(GetLayer()) {
  case (ABSLAYER):
  case (DMLAYER):
    return DPTENSORPHASE;
  case (SMLAYER):
    return NUMPHASES;
  default:
    throw;
  }
}



string DistContToLocalContStatC::GetType() const
{
  return "DistContToLocalContStatC";
}

Cost DistContToLocalContStatC::RHSCostEstimate(const Node *node) const
{
  const Contraction *cont = (Contraction*)node;
  Cost cost1 = 1;
  for (Dim dim = 0; dim < cont->InputNumDims(0); ++dim) {
    cost1 *= (*(cont->InputLen(0, dim)))[0];
  }
  Cost cost2 = 1;
  for (Dim dim = 0; dim < cont->InputNumDims(1); ++dim) {
    cost2 *= (*(cont->InputLen(1, dim)))[0];
  }
  return cost1 + cost2;
}

bool DistContToLocalContStatC::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Contraction::GetClass())
    throw;
  const Contraction *cont = (Contraction*)node;
  if (!cont->InputDataType(2).m_dist.HasNoReped())
    return false;
  return (cont->GetLayer() == m_fromLayer);
}

void DistContToLocalContStatC::Apply(Node *node) const
{
  Contraction *cont = (Contraction*)node;

  DimVec ADims = MapIndicesToDims(cont->m_contIndices,cont->m_AIndices);
  DimVec BDims = MapIndicesToDims(cont->m_contIndices,cont->m_BIndices);

  NodeConn *CConn = cont->InputConn(2);
  if (CConn->m_n->GetNodeClass() == RedistNode::GetClass())
    CConn = CConn->m_n->InputConn(0);

  const DistType &CType = ((DLANode*)(CConn->m_n))->DataType(CConn->m_num).m_dist;

  DistType AType;
  MatchDistsAndFillInWithStar(cont->m_AIndices,
			      CType, cont->m_CIndices,
			      AType);

  
  DistType BType;
  MatchDistsAndFillInWithStar(cont->m_BIndices,
			      CType, cont->m_CIndices,
			      BType);

  RedistNode *node1 = new RedistNode(AType);
  node1->AddInput(node->Input(0),node->InputConnNum(0));

  if (AType == node1->InputDataType(0).m_dist)
    throw;

  Poss *APoss = new Poss(node1, false);
  RealPSet *ASet = new RealPSet(APoss);
  node->m_poss->AddPSet(ASet,true,true);

  RedistNode *node2 = new RedistNode(BType);
  node2->AddInput(node->Input(1),node->InputConnNum(1));

  if (BType == node2->InputDataType(0).m_dist)
    throw;

  Poss *BPoss = new Poss(node2, false);
  RealPSet *BSet = new RealPSet(BPoss);
  node->m_poss->AddPSet(BSet,true,true);


  Contraction *LCont = new Contraction(m_toLayer,  cont->m_alpha, cont->m_beta, cont->m_type, 
				       cont->m_AIndices, cont->m_BIndices, cont->m_CIndices, cont->m_contIndices);
  LCont->AddInput(ASet->OutTun(0),0);
  LCont->AddInput(BSet->OutTun(0),0);
  LCont->AddInput(node->Input(2),node->InputConnNum(2));
  node->m_poss->AddNode(LCont);

  cont->RedirectChildren(LCont,0);

  node->m_poss->DeleteChildAndCleanUp(node);
}



string DistContToLocalContStatAAllReduce::GetType() const
{
  return "DistContToLocalContStatAAllReduce";
}

bool DistContToLocalContStatAAllReduce::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Contraction::GetClass())
    throw;
  const Contraction *cont = (Contraction*)node;
  if (!cont->InputDataType(0).m_dist.HasNoReped())
    return false;
  if (!cont->InputDataType(2).m_dist.HasNoReped())
    return false;
  return (cont->GetLayer() == m_fromLayer);
}

void DistContToLocalContStatAAllReduce::Apply(Node *node) const
{
  Contraction *cont = (Contraction*)node;

  DimVec BDims = MapIndicesToDims(cont->m_contIndices,cont->m_BIndices);
  DimVec CDims = MapIndicesToDims(cont->m_contIndices,cont->m_CIndices);

  NodeConn *AConn = cont->InputConn(0);
  if (AConn->m_n->GetNodeClass() == RedistNode::GetClass())
    AConn = AConn->m_n->InputConn(0);

  const DistType &AType = ((DLANode*)(AConn->m_n))->DataType(AConn->m_num).m_dist;

  DimVec sumDims;
  string sumIndices;
  string::iterator iter = cont->m_contIndices.begin();
  for(; iter != cont->m_contIndices.end(); ++iter) {
    size_t loc = cont->m_AIndices.find(*iter);
    if (loc != string::npos) {
      DimVec dims = AType.m_dists[loc].DistEntryDims();
      DimVecIter iter2 = dims.begin();
      for( ; iter2 != dims.end(); ++iter2) {
	sumDims.push_back(*iter2);
	sumIndices += *iter;
      }
    }
  }

  DistType BType;
  MatchDistsAndFillInWithStar(cont->m_BIndices,
			      AType, cont->m_AIndices,
			      BType);
  
  DistType CType;
  MatchDistsAndFillInWithStar(cont->m_CIndices,
			      AType, cont->m_AIndices, 
			      CType);

  


  RedistNode *node1 = new RedistNode(AType);
  RedistNode *node2 = new RedistNode(BType);
  RedistNode *node3 = new RedistNode(CType);
  Contraction *LCont = new Contraction(m_toLayer,  cont->m_alpha, cont->m_beta, cont->m_type, 
				       cont->m_AIndices, cont->m_BIndices, cont->m_CIndices, cont->m_contIndices);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node->Input(1),node->InputConnNum(1));
  node3->AddInput(node->Input(2),node->InputConnNum(2));
  LCont->AddInput(node1,0);
  LCont->AddInput(node2,0);
  LCont->AddInput(node3,0);
  node->m_poss->AddNode(node1);
  node->m_poss->AddNode(node2);
  node->m_poss->AddNode(node3);
  node->m_poss->AddNode(LCont);

  AllReduceNode *sum = new AllReduceNode(sumDims, sumIndices);
  sum->AddInput(LCont, 0);
  node->m_poss->AddNode(sum);

  RedistNode *node4 = new RedistNode(cont->DataType(0).m_dist);
  node4->AddInput(sum);
  node->m_poss->AddNode(node4);

  cont->RedirectChildren(node4,0);

  node->m_poss->DeleteChildAndCleanUp(node);
}




string DistContToLocalContStatASumScatter::GetType() const
{
  return "DistContToLocalContStatASumScatter";
}

bool DistContToLocalContStatASumScatter::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Contraction::GetClass())
    throw;
  const Contraction *cont = (Contraction*)node;
  if (!cont->InputDataType(0).m_dist.HasNoReped())
    return false;
  return (cont->GetLayer() == m_fromLayer);
}

Cost DistContToLocalContStatASumScatter::RHSCostEstimate(const Node *node) const
{
  const Contraction *cont = (Contraction*)node;
  Cost cost1 = 1;
  for (Dim dim = 0; dim < cont->InputNumDims(1); ++dim) {
    cost1 *= (*(cont->InputLen(1, dim)))[0];
  }
  Cost cost2 = 1;
  for (Dim dim = 0; dim < cont->InputNumDims(2); ++dim) {
    cost2 *= (*(cont->InputLen(2, dim)))[0];
  }
  return cost1 + cost2;
}

void DistContToLocalContStatASumScatter::Apply(Node *node) const
{
  Contraction *cont = (Contraction*)node;

  DimVec BDims = MapIndicesToDims(cont->m_contIndices,cont->m_BIndices);
  DimVec CDims = MapIndicesToDims(cont->m_contIndices,cont->m_CIndices);


  NodeConn *AConn = cont->InputConn(0);
  if (AConn->m_n->GetNodeClass() == RedistNode::GetClass())
    AConn = AConn->m_n->InputConn(0);

  const DistType &AType = ((DLANode*)(AConn->m_n))->DataType(AConn->m_num).m_dist;

  EntryList sumDims;
  string::iterator iter = cont->m_contIndices.begin();
  for(; iter != cont->m_contIndices.end(); ++iter) {
    size_t loc = cont->m_AIndices.find(*iter);
    if (loc != string::npos) {
      DistEntry entry = AType.m_dists[loc];
      sumDims.push_back(entry);
    }
  }

  DistType BType;
  MatchDistsAndFillInWithStar(cont->m_BIndices,
			      AType, cont->m_AIndices,
			      BType);

  RedistNode *node1 = NULL;
  RealPSet *BSet = NULL;

  if (BType != node->InputDataType(1).m_dist) {
    node1 = new RedistNode(BType);
    node1->AddInput(node->Input(1),node->InputConnNum(1));

    Poss *BPoss = new Poss(node1, false);
    BSet = new RealPSet(BPoss);
    node->m_poss->AddPSet(BSet,true,true);
  }

  
  DistType CType;
  const DistType &CDestType = cont->InputDataType(2).m_dist;
  MatchDistsAndFillInWithStar(cont->m_CIndices,
			      AType, cont->m_AIndices, 
			      CType);

  TempVarNode *temp = new TempVarNode(CType, sumDims);

  Contraction *LCont = new Contraction(m_toLayer,  cont->m_alpha, COEFVALZERO, cont->m_type, 
				       cont->m_AIndices, cont->m_BIndices, cont->m_CIndices+cont->m_contIndices, cont->m_contIndices);
  temp->AddInput(node->Input(2),node->InputConnNum(2));
  LCont->AddInput(node->Input(0),node->InputConnNum(0));
  if (BSet)
    LCont->AddInput(BSet->OutTun(0),0);
  else
    LCont->AddInput(node->Input(1),node->InputConnNum(1));
  LCont->AddInput(temp,0);

  node->m_poss->AddNode(temp);
  node->m_poss->AddNode(LCont);

  bool notFinalType = false;
  iter = cont->m_CIndices.begin();
  EntryListIter iter2 = sumDims.begin();
  for(int i = 0; iter != cont->m_CIndices.end(); ++iter, ++i) {
    size_t loc = cont->m_AIndices.find(*iter);
    if (loc != string::npos) {
      CType.m_dists[i] = AType.m_dists[loc]; 
      if (CType.m_dists[i] != CDestType.m_dists[i])
        notFinalType = true;
    }
    else {
      if (iter2 == sumDims.end())
	throw;
      CType.m_dists[i] = *iter2;
      ++iter2;
    }
  }


  if (!notFinalType) {
    SumScatterUpdateNode *sum = new SumScatterUpdateNode(cont->m_beta, sumDims);
    sum->AddInput(LCont, 0);
    sum->AddInput(node->Input(2),node->InputConnNum(2));
    
    Poss *sumPoss = new Poss(sum, false);
    RealPSet *sumSet = new RealPSet(sumPoss);
    node->m_poss->AddPSet(sumSet,true,true);
    cont->RedirectChildren(sumSet->OutTun(0),0);
  }
  else {
    TempVarNode *temp2;
    if (cont->m_beta != COEFZERO)
	temp2 = new TempVarNode(CType,node->GetInputName(2).m_name+"_temp");
    else
      temp2 = new TempVarNode(CType);
    temp2->AddInput(node->Input(2),node->InputConnNum(2));
    node->m_poss->AddNode(temp2);

    SumScatterUpdateNode *sum = new SumScatterUpdateNode(COEFZERO, sumDims);
    sum->AddInput(LCont, 0);
    sum->AddInput(temp2, 0);
    Poss *sumPoss = new Poss(sum, false);
    RealPSet *sumSet = new RealPSet(sumPoss);
    node->m_poss->AddPSet(sumSet,true,true);
    
    RedistNode *finalRedist = new RedistNode(CDestType);
    finalRedist->AddInput(sumSet->OutTun(0),0);
    Poss *redistPoss = new Poss(finalRedist, false);
    RealPSet *redistSet = new RealPSet(redistPoss);
    node->m_poss->AddPSet(redistSet,true,true);
    
    if (cont->m_beta != COEFZERO) {
      YxpBy *yxpby = new YxpBy(SMLAYER, cont->m_beta);
      yxpby->AddInput(redistSet->OutTun(0), 0);
      yxpby->AddInput(node->Input(2),node->InputConnNum(2));
      node->m_poss->AddNode(yxpby);
      cont->RedirectChildren(yxpby,0);
    }
    else {
      cont->RedirectChildren(redistSet->OutTun(0),0);
    }
  }


  node->m_poss->DeleteChildAndCleanUp(node);
}





string DistContToLocalContStatBAllReduce::GetType() const
{
  return "DistContToLocalContStatBAllReduce";
}

bool DistContToLocalContStatBAllReduce::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Contraction::GetClass())
    throw;
  const Contraction *cont = (Contraction*)node;
  if (!cont->InputDataType(1).m_dist.HasNoReped())
    return false;
  if (!cont->InputDataType(2).m_dist.HasNoReped())
    return false;
  return (cont->GetLayer() == m_fromLayer);
}

void DistContToLocalContStatBAllReduce::Apply(Node *node) const
{
  Contraction *cont = (Contraction*)node;

  DimVec ADims = MapIndicesToDims(cont->m_contIndices,cont->m_AIndices);
  DimVec CDims = MapIndicesToDims(cont->m_contIndices,cont->m_CIndices);


  NodeConn *BConn = cont->InputConn(1);
  if (BConn->m_n->GetNodeClass() == RedistNode::GetClass())
    BConn = BConn->m_n->InputConn(0);

  const DistType &BType = ((DLANode*)(BConn->m_n))->DataType(BConn->m_num).m_dist;

  DimVec sumDims;
  string sumIndices;
  string::iterator iter = cont->m_contIndices.begin();
  for(; iter != cont->m_contIndices.end(); ++iter) {
    size_t loc = cont->m_BIndices.find(*iter);
    if (loc != string::npos) {
      DimVec dims = BType.m_dists[loc].DistEntryDims();
      DimVecIter iter2 = dims.begin();
      for(; iter2 != dims.end(); ++iter2) {
	sumDims.push_back(*iter2);
	sumIndices += *iter;
      }
    }
  }

  DistType AType;
  MatchDistsAndFillInWithStar(cont->m_AIndices,
			      BType, cont->m_BIndices,
			      AType);
  
  DistType CType;
  MatchDistsAndFillInWithStar(cont->m_CIndices,
			      BType, cont->m_BIndices, 
			      CType);

  


  RedistNode *node1 = new RedistNode(AType);
  RedistNode *node2 = new RedistNode(BType);
  RedistNode *node3 = new RedistNode(CType);
  Contraction *LCont = new Contraction(m_toLayer,  cont->m_alpha, cont->m_beta, cont->m_type, 
				       cont->m_AIndices, cont->m_BIndices, cont->m_CIndices, cont->m_contIndices);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node->Input(1),node->InputConnNum(1));
  node3->AddInput(node->Input(2),node->InputConnNum(2));
  LCont->AddInput(node1,0);
  LCont->AddInput(node2,0);
  LCont->AddInput(node3,0);
  node->m_poss->AddNode(node1);
  node->m_poss->AddNode(node2);
  node->m_poss->AddNode(node3);
  node->m_poss->AddNode(LCont);

  AllReduceNode *sum = new AllReduceNode(sumDims, sumIndices);
  sum->AddInput(LCont, 0);
  node->m_poss->AddNode(sum);

  RedistNode *node4 = new RedistNode(cont->DataType(0).m_dist);
  node4->AddInput(sum);
  node->m_poss->AddNode(node4);

  cont->RedirectChildren(node4,0);

  node->m_poss->DeleteChildAndCleanUp(node);
}



Cost DistContToLocalContStatBSumScatter::RHSCostEstimate(const Node *node) const
{
  const Contraction *cont = (Contraction*)node;
  Cost cost1 = 1;
  for (Dim dim = 0; dim < cont->InputNumDims(0); ++dim) {
    cost1 *= (*(cont->InputLen(0, dim)))[0];
  }
  Cost cost2 = 1;
  for (Dim dim = 0; dim < cont->InputNumDims(2); ++dim) {
    cost2 *= (*(cont->InputLen(2, dim)))[0];
  }
  return cost1 + cost2;
}


string DistContToLocalContStatBSumScatter::GetType() const
{
  return "DistContToLocalContStatBSumScatter";
}

bool DistContToLocalContStatBSumScatter::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Contraction::GetClass())
    throw;
  const Contraction *cont = (Contraction*)node;
  if (!cont->InputDataType(1).m_dist.HasNoReped())
    return false;
  return (cont->GetLayer() == m_fromLayer);
}

void DistContToLocalContStatBSumScatter::Apply(Node *node) const
{
  Contraction *cont = (Contraction*)node;

  DimVec ADims = MapIndicesToDims(cont->m_contIndices,cont->m_AIndices);
  DimVec CDims = MapIndicesToDims(cont->m_contIndices,cont->m_CIndices);


  NodeConn *BConn = cont->InputConn(1);
  if (BConn->m_n->GetNodeClass() == RedistNode::GetClass())
    BConn = BConn->m_n->InputConn(0);

  const DistType &BType = ((DLANode*)(BConn->m_n))->DataType(BConn->m_num).m_dist;

  EntryList sumDims;
  string::iterator iter = cont->m_contIndices.begin();
  for(; iter != cont->m_contIndices.end(); ++iter) {
    size_t loc = cont->m_BIndices.find(*iter);
    if (loc != string::npos) {
      DistEntry entry = BType.m_dists[loc];
      sumDims.push_back(entry);
    }
  }

  DistType AType;
  MatchDistsAndFillInWithStar(cont->m_AIndices,
			      BType, cont->m_BIndices,
			      AType);

  RedistNode *node1 = NULL;
  RealPSet *ASet = NULL;

  if (AType != node->InputDataType(0).m_dist) {
    node1 = new RedistNode(AType);
    node1->AddInput(node->Input(0),node->InputConnNum(0));

    Poss *APoss = new Poss(node1, false);
    ASet = new RealPSet(APoss);
    node->m_poss->AddPSet(ASet,true,true);
  }

  //  cout << "AType " << AType.PrettyStr() << endl;

  
  DistType CType;
  const DistType &CDestType = cont->InputDataType(2).m_dist;
  MatchDistsAndFillInWithStar(cont->m_CIndices,
			      BType, cont->m_BIndices, 
			      CType);

  TempVarNode *temp = new TempVarNode(CType, sumDims);

  Contraction *LCont = new Contraction(m_toLayer,  cont->m_alpha, COEFVALZERO, cont->m_type, 
				       cont->m_AIndices, cont->m_BIndices, cont->m_CIndices+cont->m_contIndices, cont->m_contIndices);
  temp->AddInput(node->Input(2),node->InputConnNum(2));
  if (ASet)
    LCont->AddInput(ASet->OutTun(0),0);
  else
    LCont->AddInput(node->Input(0),node->InputConnNum(0));
  LCont->AddInput(node->Input(1),node->InputConnNum(1));
  LCont->AddInput(temp,0);

  node->m_poss->AddNode(temp);
  node->m_poss->AddNode(LCont);

  bool notFinalType = false;
  iter = cont->m_CIndices.begin();
  EntryListIter iter2 = sumDims.begin();
  for(int i = 0; iter != cont->m_CIndices.end(); ++iter, ++i) {
    size_t loc = cont->m_BIndices.find(*iter);
    if (loc != string::npos) {
      CType.m_dists[i] = BType.m_dists[loc]; 
      if (CType.m_dists[i] != CDestType.m_dists[i])
        notFinalType = true;
    }
    else {
      if (iter2 == sumDims.end())
	throw;
      CType.m_dists[i] = *iter2;
      ++iter2;
    }
  }

  if (!notFinalType) {
    SumScatterUpdateNode *sum = new SumScatterUpdateNode(cont->m_beta, sumDims);
    sum->AddInput(LCont, 0);
    sum->AddInput(node->Input(2),node->InputConnNum(2));
    
    Poss *sumPoss = new Poss(sum, false);
    RealPSet *sumSet = new RealPSet(sumPoss);
    node->m_poss->AddPSet(sumSet,true,true);
    cont->RedirectChildren(sumSet->OutTun(0),0);
  }
  else {
    TempVarNode *temp2;
    if (cont->m_beta != COEFZERO)
      temp2 = new TempVarNode(CType, node->GetInputName(2).m_name+"_temp");
    else
      temp2 = new TempVarNode(CType);
    temp2->AddInput(node->Input(2),node->InputConnNum(2));
    node->m_poss->AddNode(temp2);

    SumScatterUpdateNode *sum = new SumScatterUpdateNode(COEFZERO, sumDims);
    sum->AddInput(LCont, 0);
    sum->AddInput(temp2, 0);
    Poss *sumPoss = new Poss(sum, false);
    RealPSet *sumSet = new RealPSet(sumPoss);
    node->m_poss->AddPSet(sumSet,true,true);

    RedistNode *finalRedist = new RedistNode(CDestType);
    finalRedist->AddInput(sumSet->OutTun(0),0);
    Poss *redistPoss = new Poss(finalRedist, false);
    RealPSet *redistSet = new RealPSet(redistPoss);
    node->m_poss->AddPSet(redistSet,true,true);
    
    if (cont->m_beta != COEFZERO) {
      YxpBy *yxpby = new YxpBy(SMLAYER, cont->m_beta);
      yxpby->AddInput(redistSet->OutTun(0), 0);
      yxpby->AddInput(node->Input(2),node->InputConnNum(2));
      node->m_poss->AddNode(yxpby);
      cont->RedirectChildren(yxpby,0);
    }
    else {
      cont->RedirectChildren(redistSet->OutTun(0),0);
    }
  }

  
  node->m_poss->DeleteChildAndCleanUp(node);
}
  
void MatchDistsAndFillInWithStar(string indices, 
				 const DistType &matchingDists, string matchingIndices,
				 DistType &final)
{
  final.PrepForNumDims(indices.length());
  for(unsigned int i = 0; i < indices.length(); ++i) {
    char index = indices[i];
    size_t loc = matchingIndices.find(index);
    if (loc == string::npos) {
      final.m_dists[i].SetToStar();
    }
    else {
      final.m_dists[i] = matchingDists.m_dists[loc];
    }
  }
}



/*
int DistContToLocalContStatC::CanApply(const Node *node, void **cache) const
{
  if (node->GetNodeClass() != Contraction::GetClass())
    throw;
  const Contraction *cont = (Contraction*)node;
  
  NodeConn *CConn = cont->InputConn(2);
  if (CConn->m_n->GetNodeClass() == RedistNode::GetClass())
    CConn = CConn->m_n->InputConn(0);
  const DistType &CType = ((DLANode*)(CConn->m_n))->GetDistType(CConn->m_num);
  if (!CType.m_numDims)
    return false;
  
  Dim numContDims = cont->m_contIndices.length();
  
  DimVec ADims = MapIndicesToDims(cont->m_contIndices,cont->m_AIndices);
  DimVec BDims = MapIndicesToDims(cont->m_contIndices,cont->m_BIndices);

  NodeConn *AConn = cont->InputConn(0);
  if (AConn->m_n->GetNodeClass() == RedistNode::GetClass())
    AConn = AConn->m_n->InputConn(0);
  NodeConn *BConn = cont->InputConn(1);
  if (BConn->m_n->GetNodeClass() == RedistNode::GetClass())
    BConn = BConn->m_n->InputConn(0);



  const DistType &AType = ((DLANode*)(AConn->m_n))->GetDistType(AConn->m_num);
  const DistType &BType = ((DLANode*)(BConn->m_n))->GetDistType(BConn->m_num);


  DimVec *dists = new DimVec[numContDims];
  
  DimSet usedDims = CType.UsedGridDims();

  ContType *distOptions = new ContType;

  RecursivelyFindDistributions(dists, 0, AType, ADims, BType, BDims, usedDims, CType.m_numDims+numContDims, distOptions);

  *cache = distOptions;
  
  delete [] dists;

  return distOptions->size();
}

void RecursivelyFindDistributions(DimVec *dists, Dim thisDim, 
				  const DistType &AType, const DimVec &ADims,
				  const DistType &BType, const DimVec &BDims,
				  DimSet &usedDims, Dim numDims,
				  ContType *distOptions)
{

//     dists ->  a numContDims-length c-style array of DimVecs;
//               the kth DimVec holds the dimensions for
// 	      the k-th index's distribution
//     thisDim -> this call is the thisDim^{th} index's call
//     {A,B}Type -> DistType for the {A,B} tensor
//     {A,B}Dims -> Map of contraction indices to dimensions of {A,B} (i.e.,
//                  where those indices are in the tensor}
// 		 So the ADim[thisDim] distribution of AType is the
// 		 distribution of A for the current contraction index
//     usedDims -> List of processes grid dimensions that have been 
//                 used for distribution up to this point in the recursion
//     distOptions -> Vector of DistTypes 
  
  if (thisDim == ADims.size()) {
    DimSet unUsedDims;
    DimSetIter iter = usedDims.begin();
    for(Dim dim = 0; dim < numDims; ++dim) {
      if (iter != usedDims.end()) {
	if (*iter > dim) {
	  unUsedDims.insert(dim);
	}
	else if (*iter == dim)
	  ++iter;
	else
	  throw;
      }
      else
	unUsedDims.insert(dim);
    }
    if (iter != usedDims.end())
      throw;
    unsigned int *distEntries = new unsigned int [ADims.size()];
    for (Dim dim = 0; dim < ADims.size(); ++dim)
      distEntries[dim] = 0;
    AddUnusedDimsForDistType(dists, distEntries, ADims.size(), unUsedDims, distOptions);
    delete [] distEntries;
    return;
  }
  //Fist, call recursively with * for this dim
  RecursivelyFindDistributions(dists, thisDim+1, 
			       AType, ADims,
			       BType, BDims,
			       usedDims, numDims,
			       distOptions);
  // Now, call recursively with A's Type (if it's not already used)
  DimVec ADists = DistType::DistEntryDims(AType.m_dists[ADims[thisDim]]);
  DimVecIter AIter = ADists.begin();
  DimSet usedDimsTemp = usedDims;
  for(; AIter != ADists.end(); ++AIter) {
    Dim dim = *AIter;
    //check that this isn't a distribution we've already used
    if (usedDims.find(dim) != usedDims.end())
      break;
    usedDimsTemp.insert(dim);
    dists[thisDim].push_back(dim);
    RecursivelyFindDistributions(dists, thisDim+1,
				 AType, ADims,
				 BType, BDims,
				 usedDimsTemp, numDims,
				 distOptions);
  }

  dists[thisDim].clear();
  AIter = ADists.begin();
  DimVec BDists = DistType::DistEntryDims(BType.m_dists[BDims[thisDim]]);
  DimVecIter BIter = BDists.begin();
  usedDimsTemp = usedDims;
  bool stillMatchingA = true;
  for(; BIter != BDists.end(); ++BIter, ++AIter) {
    Dim dim = *BIter;
    if (stillMatchingA) {
      if (dim == *AIter) {
	usedDimsTemp.insert(dim);
	dists[thisDim].push_back(dim);
	continue;
      }
      else
	stillMatchingA = false;
    }
    //check that this isn't a distribution we've already used
    if (usedDims.find(dim) != usedDims.end())
      break;
    usedDimsTemp.insert(dim);
    dists[thisDim].push_back(dim);
    RecursivelyFindDistributions(dists, thisDim+1,
				 AType, ADims,
				 BType, BDims,
				 usedDimsTemp, numDims,
				 distOptions);
  }
}

void AddUnusedDimsForDistType(DimVec *dists,  unsigned int *distEntries,
			      Dim numIndices,
			       DimSet &unUsedDims,
			       ContType *distOptions)
{
  unsigned int *entries = new unsigned int[numIndices];
  for (Dim dim = 0; dim < numIndices; ++dim) {
    if (distEntries[dim] == 0 && !dists[dim].empty()) {
      distEntries[dim] = DistType::DimsToDistEntry(dists[dim]);
    }
    entries[dim] = distEntries[dim];
  }
  distOptions->push_back(entries);
  if (unUsedDims.empty()) {
    return;
  }
  DimSetIter iter = unUsedDims.begin();
  for(; iter != unUsedDims.end(); ++iter) {
    Dim unused = *iter;
    DimSet tempUnUsedDims = unUsedDims;
    tempUnUsedDims.erase(unused);
    for (Dim dim = 0; dim < numIndices; ++dim) {
      unsigned int temp = distEntries[dim];
      distEntries[dim] = 0;
      dists[dim].push_back(unused);
      AddUnusedDimsForDistType(dists, distEntries,
			       numIndices,
			       tempUnUsedDims,
			       distOptions);
      dists[dim].pop_back();
      distEntries[dim] = temp;
    }
  }
}

void DistContToLocalContStatC::Apply(int num, Node *node, void **cache) const
{
  ContType *types = (ContType*)(*cache);
  unsigned int *entries = (*types)[num];

  Contraction *cont = (Contraction*)node;
  
  DimVec ADims = MapIndicesToDims(cont->m_contIndices,cont->GetInputName(0).m_indices);
  DimVec BDims = MapIndicesToDims(cont->m_contIndices,cont->GetInputName(1).m_indices);

  NodeConn *CConn = cont->InputConn(2);
  if (CConn->m_n->GetNodeClass() == RedistNode::GetClass())
    CConn = CConn->m_n->InputConn(0);

    const DistType &CType = ((DLANode*)(CConn->m_n))->DataType(CConn->m_num).m_dist;

  DistType AType;
  MatchDistsAndFillIn(cont->GetInputName(0).m_indices,
		      CType, ((DLANode*)(CConn->m_n))->GetName(CConn->m_num).m_indices,
		      entries, ADims,
		      AType);
  
  DistType BType;
  MatchDistsAndFillIn(cont->GetInputName(1).m_indices,
		      CType, ((DLANode*)(CConn->m_n))->GetName(CConn->m_num).m_indices,
		      entries, BDims,
		      BType);

  RedistNode *node1 = new RedistNode(AType);
  node1->AddInput(node->Input(0),node->InputConnNum(0));

  Poss *APoss = new Poss(node1, false);
  RealPSet *ASet = new RealPSet(APoss);
  node->m_poss->AddPSet(ASet,true,true);

  RedistNode *node2 = new RedistNode(BType);
  node2->AddInput(node->Input(1),node->InputConnNum(1));

  Poss *BPoss = new Poss(node2, false);
  RealPSet *BSet = new RealPSet(BPoss);
  node->m_poss->AddPSet(BSet,true,true);

  RedistNode *node3 = new RedistNode(CType);
  node3->AddInput(node->Input(2),node->InputConnNum(2));

  Poss *CPoss = new Poss(node3, false);
  RealPSet *CSet = new RealPSet(CPoss);
  node->m_poss->AddPSet(CSet,true,true);

  Contraction *LCont = new Contraction(SMLAYER,  cont->m_alpha, cont->m_beta, cont->m_type, cont->m_indices);
  LCont->AddInput(ASet->OutTun(0),0);
  LCont->AddInput(BSet->OutTun(0),0);
  LCont->AddInput(CSet->OutTun(0),0);

  bool sum = false;
  for (unsigned int i = 0; i < cont->m_indices.length() && !sum; ++i) {
    if (entries[i] != 0)
      sum = true;
  }

  DLANode *node4;
  if (sum)
    cout << "need different refinement code!\n";
  //  else {
    node4 = new RedistNode(cont->InputDataType(2).m_dist);
    node4->AddInput(LCont, 0);

    Poss *CAfterPoss = new Poss(node4, false);
    RealPSet *CAfterSet = new RealPSet(CAfterPoss);
    node->m_poss->AddPSet(CAfterSetSet,true,true);
    //  }
    
  node->m_poss->AddNode(LCont);

  node->RedirectChildren(CAfterSet->OutTun(0),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

void MatchDistsAndFillIn(string indices, 
			 const DistType &matchingDists, string matchingIndices,
			 unsigned int *fillInDists, const DimVec &fillInDims,
			 DistType &final)
{
  unsigned int tot = 0;
  final.PrepForNumDims(indices.length());
  string::iterator iter = indices.begin();
  for(; iter != indices.end(); ++iter) {
    char index = *iter;
    size_t loc = matchingIndices.find(index);
    if (loc != string::npos) {
      final.m_dists[loc] = matchingDists.m_dists[loc];
      ++tot;
    }
  }
  if ((tot + fillInDims.size()) != indices.length()) {
    cout << "tot = " << tot << endl;
    cout << "fillInDims.size() = " << fillInDims.size() << endl;
    cout << "indices.length() = " << indices.length() << endl;
    throw;
  }
  DimVecConstIter iter2 = fillInDims.begin();
  for(Dim dimNum = 0; iter2 != fillInDims.end(); ++iter2, ++dimNum) {
    Dim dim = *iter2;
    final.m_dists[dim] = fillInDists[dimNum];
  }
}

void DistContToLocalContStatC::CleanCache(void **cache) const
{
  ContType *types = (ContType*)(*cache);
  ContTypeIter iter = types->begin();
  for(; iter != types->end(); ++iter)
    delete [] *iter;
  delete types;
}
*/
/*
Cost DistContToLocalContStatC::RHSCostEstimate(const Node *node) const
{
  Cost cost = 0;
  const Contraction *cont = (Contraction*)node;
  DimVec dims = MapIndicesToDims(cont->m_indices,cont->GetInputName(0).m_indices);
  const Sizes *sizes = cont->InputLocalLen(2,0);
  unsigned int totNumIters = sizes->NumSizes();
  for(unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
    Cost temp = 1;
    Dim numDims = cont->InputNumDims(2);
    for (Dim dim = 1; dim < numDims; ++dim) {
      temp *= (*(cont->InputLocalLen(2,dim)))[iteration];
    }
    DimVecConstIter iter = dims.begin();
    for(; iter != dims.end(); ++iter) {
      temp *= (*(cont->InputLocalLen(0,*iter)))[iteration];
    }
    cost += temp;
  }
  return cost * 2;
}
*/


ContractionLoopExp::ContractionLoopExp(Layer fromLayer, Layer toLayer, int dim)
  : 
  m_fromLayer(fromLayer), 
  m_toLayer(toLayer), 
  m_dim(dim) 
{
}

string ContractionLoopExp::GetType() const
{
  string str = "Contraction Loop Exp " 
    + LayerNumToStr(m_fromLayer)
    + " + " 
    + LayerNumToStr(m_toLayer)
    + " on dim " + std::to_string(m_dim);
  return str;
}

bool ContractionLoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Contraction::GetClass()) {
    const Contraction *cont = (Contraction*)node;
    if (cont->GetLayer() == m_fromLayer) {
      const string &cIndices = cont->m_CIndices;
      if (cIndices.length() > m_dim) {
	return !(*(cont->InputLen(2,m_dim)) <= TensorBS.GetSize());
      }
      else if (cont->m_contIndices.length() > (m_dim - cIndices.length())) {
	char contIndex = cont->m_contIndices[m_dim - cIndices.length()];
	size_t index = cont->m_AIndices.find(contIndex);
	if (index == string::npos)
	  throw;
	return !(*(cont->InputLen(0,index)) <= TensorBS.GetSize());
      }
    }
  }
  return false;
}

void ContractionLoopExp::Apply(Node *node) const
{
  Contraction *cont = (Contraction*)node;
  const string &aIndices = cont->m_AIndices;
  const string &bIndices = cont->m_BIndices;
  const string &cIndices = cont->m_CIndices;
  const bool isCIndex = cIndices.length() > m_dim;
  const char splitIndex = (isCIndex 
			   ? cIndices[m_dim] 
			   : cont->m_contIndices[m_dim-cIndices.length()]);
  const size_t cDim = m_dim;

  const size_t aDim = aIndices.find(splitIndex);
  const bool isAIndex = aDim != string::npos;
  const size_t bDim = bIndices.find(splitIndex);
  const bool isBIndex = bDim != string::npos;
  

  NodeConn *connA, *connB, *connC;
  connA = cont->m_inputs[0];
  connB = cont->m_inputs[1];
  connC = cont->m_inputs[2];

  LoopTunnel *ATun;
  ConnNum ATunNum;
  if (isAIndex) {
    ATun = new SplitSingleIter(aDim, POSSTUNIN, isCIndex ? false : true);
    ATunNum = 1;
  }
  else {
    ATun = new LoopTunnel(POSSTUNIN);
    ATunNum = 0;
  }
  ATun->AddInput(connA->m_n, connA->m_num);
  ATun->SetAllStats(FULLUP);
  if (isCIndex)
    ATun->SetIndepIters();


  LoopTunnel *BTun;
  ConnNum BTunNum;
  if (isBIndex) {
    BTun = new SplitSingleIter(bDim, POSSTUNIN, false);
    BTunNum = 1;
  }
  else {
    BTun = new LoopTunnel(POSSTUNIN);
    BTunNum = 0;
  }
  BTun->AddInput(connB->m_n, connB->m_num);
  BTun->SetAllStats(FULLUP);
  if (isCIndex)
    BTun->SetIndepIters();


  LoopTunnel *CTun;
  ConnNum CTunNum;
  if (isCIndex) {
    CTun = new SplitSingleIter(cDim, POSSTUNIN, true);
    CTunNum = 1;
    CTun->AddInput(connC->m_n, connC->m_num);
  }
  else {
    CTun = new LoopTunnel(POSSTUNIN);
    CTunNum = 0;
    if (cont->m_beta == COEFONE)
      CTun->AddInput(connC->m_n, connC->m_num);
    else {
      ScaleNode *scale = new ScaleNode(m_toLayer, cont->m_beta);
      cont->m_poss->AddNode(scale);
      scale->AddInput(connC->m_n, connC->m_num);
      CTun->AddInput(scale,0);
    }
  }
  CTun->SetAllStats(FULLUP);
  if (isCIndex)
    CTun->SetIndepIters();
  
  Contraction *newCont = new Contraction(m_toLayer, cont->m_alpha, 
					 isCIndex ? cont->m_beta : COEFZERO, 
					 cont->m_type,
					 cont->m_AIndices,
					 cont->m_BIndices,
					 cont->m_CIndices,
					 cont->m_contIndices);
  newCont->AddInput(ATun, ATunNum);
  newCont->AddInput(BTun, BTunNum);
  newCont->AddInput(CTun, CTunNum);
  
  LoopTunnel *AOut;
  if (isAIndex)
    AOut = ((SplitSingleIter*)ATun)->CreateMatchingCombine(0);
  else {
    AOut = new LoopTunnel(POSSTUNOUT);
    AOut->AddInput(ATun, 0);
    AOut->AddInput(ATun, 1);
    AOut->CopyTunnelInfo(ATun);
  }


  LoopTunnel *BOut;
  if (isBIndex)
    BOut = ((SplitSingleIter*)BTun)->CreateMatchingCombine(0);
  else {
    BOut = new LoopTunnel(POSSTUNOUT);
    BOut->AddInput(BTun, 0);
    BOut->AddInput(BTun, 1);
    BOut->CopyTunnelInfo(BTun);
  }


  LoopTunnel *COut;
  if (isCIndex) {
    COut = ((SplitSingleIter*)CTun)->CreateMatchingCombine(1,
							   1, newCont, 0);
  }
  else {
    COut = new LoopTunnel(POSSTUNOUT);
    COut->AddInput(newCont, 0);
    COut->AddInput(CTun, 1);
    COut->CopyTunnelInfo(CTun);
  }
  
					
  
  Poss *loopPoss = new Poss(3, AOut, BOut, COut);
  RealLoop *loop = new RealLoop(TENSORLOOP, loopPoss, TensorBS);
  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(2),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#endif


