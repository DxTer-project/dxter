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
#include "tensorPermute.h"

typedef vector<unsigned int *> ContType;
typedef ContType::iterator ContTypeIter;


void MatchDistsAndFillInWithStar(string indices, 
				 const DistType &matchingDists, string matchingIndices,
				 DistType &final, DimVec &alignModes, DimVec &alignModesSrc);
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
  m_contIndices(contIndices),
  m_needsPacking(true)
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
    + LayerNumToStr(GetLayer()) + (m_needsPacking ? "withPacking" : "withoutPacking");
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
  m_needsPacking = cont->m_needsPacking;
}

void Contraction::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_alpha);
  WRITE(m_beta);
  WRITE(m_type);
  WRITE(m_needsPacking);
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
  READ(m_needsPacking);
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

    if(InputNumDims(0) != m_AIndices.size()) {
      cout << InputNumDims(0) << endl;
      cout << m_AIndices << endl;
      throw;
    }

    if(InputNumDims(1) != m_BIndices.size())
      throw;

    if(InputNumDims(2) != m_CIndices.size()) {
      cout << "num " << InputNumDims(2) << endl;
      cout << "dist " << InputDataType(2).GetEffectiveDist().PrettyStr() << endl;
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

    if (m_layer == ABSLAYER || m_layer == DM1LAYER || m_layer == DM2LAYER) {
      if (!m_needsPacking)
	throw;
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
      if (!m_needsPacking) {
	StringIter aIter = m_AIndices.begin();
	StringIter bIter = m_BIndices.begin();
	StringIter cIter = m_CIndices.begin();
	while (aIter != m_AIndices.end() 
	       && cIter != m_CIndices.end()
	       && *aIter == *cIter) 
	  {
	    if (m_contIndices.find(*aIter) != string::npos)
	      break;
	    ++aIter;
	    ++cIter;
	  }

	while (aIter != m_AIndices.end() 
	       && bIter != m_BIndices.end()
	       && *aIter == *bIter) 
	  {
	    ++aIter;
	    ++bIter;
	  }
	while (bIter != m_BIndices.end()
	       && cIter != m_CIndices.end()
	       && *bIter == *cIter) 
	  {
	    ++bIter;
	    ++cIter;
	  }
	if (cIter != m_CIndices.end()) {
	  StringIter contIter = m_contIndices.begin();
	  while (cIter != m_CIndices.end()
		 && contIter != m_contIndices.end()
		 && *cIter == *contIter) 
	    {
	      ++cIter;
	      ++contIter;
	    }
	  if (contIter != m_contIndices.end()) {
	    cout << m_AIndices << endl;
	    cout << m_BIndices << endl;
	    cout << m_CIndices << endl;
	    cout << m_contIndices << endl;
	    throw;
	  }
	}
	if (aIter != m_AIndices.end()
	    || bIter != m_BIndices.end()
	    || cIter != m_CIndices.end())
	  throw;
      }
      for(unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
	Cost temp = 1;
	for (Dim dim = 0; dim < numDims; ++dim) {
	  temp *= (*InputLocalLen(2,dim))[iteration];
	}
	if (m_needsPacking) {
	  m_cost += 2 * (PSIW + PSIR) * temp;
	}
	DimVecConstIter iter = dims.begin();
	for(; iter != dims.end(); ++iter) {
	  temp *= (*InputLocalLen(0,*iter))[iteration];
	}
	m_cost += temp * 2;
      }
      if (m_needsPacking) {
	m_cost += (PSIW + PSIR) * ((DLANode*)Input(0))->TotalNumberOfLocalElements(InputConnNum(0));
	m_cost += (PSIW + PSIR) * ((DLANode*)Input(1))->TotalNumberOfLocalElements(InputConnNum(1));
	//Input(2) handled above
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

  const DistType &AType = InputDataType(0).GetEffectiveDist();
  const DistType &BType = InputDataType(1).GetEffectiveDist();
  const DistType &CType = InputDataType(2).GetEffectiveDist();
    
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
	    *InputLocalLen(2,dim) != *InputLocalLen(1,bFind)) 
	  {
	    
	    InputLocalLen(2,dim)->Print();
	    cout << " vs. \n";
	    InputLocalLen(1,bFind)->Print();
	    cout << "C and B sizes don't align\n";

	    InputLen(2,dim)->Print();
	    cout << " vs. \n";
	    InputLen(1,bFind)->Print();

	    cout << dim << endl;
	    cout << InputDataType(2).GetEffectiveDist().PrettyStr() << endl;
	    cout << bFind << endl;
	    cout << InputDataType(1).GetEffectiveDist().PrettyStr() << endl;


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
        cout << GetInputNameStr(0) << endl;
        cout << m_AIndices << endl;
        cout << Input(0)->GetType() << endl;
        cout << GetInputNameStr(2) << endl;
        cout << m_CIndices << endl;
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
       << IndexArrayVarName(m_AIndices) 
       << (m_needsPacking ? ", true,\n" : ", false,\n");
  out.Indent(1);
  *out << in1.str()
       << ".LockedTensor(), "
       << IndexArrayVarName(m_BIndices)        
       << (m_needsPacking ? ", true,\n" : ", false,\n");
  out.Indent(1);
  out << m_beta;
  *out << ", " << in2.str() << ".Tensor(), "
       << IndexArrayVarName(m_CIndices) 
       << (m_needsPacking ? ", true);\n" : ", false);\n");

  
}


Phase Contraction::MaxPhase() const
{
  switch(GetLayer()) {
  case (ABSLAYER):
  case (DM1LAYER):
  case (DM2LAYER):
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
  if (!cont->InputDataType(2).GetDist().HasNoReped())
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

  const DistType &CType = ((DLANode*)(CConn->m_n))->DataType(CConn->m_num).GetEffectiveDist();

  DimVec alignModes, alignModesSrc;

  DistType AType;
  MatchDistsAndFillInWithStar(cont->m_AIndices,
			      CType, cont->m_CIndices,
			      AType, alignModes, alignModesSrc);

  
  RedistNode *node1 = new RedistNode(AType, cont->GetInputNameStr(2), alignModes, alignModesSrc);
  node1->AddInput(node->Input(0),node->InputConnNum(0));

  if (AType == node1->InputDataType(0).GetEffectiveDist())
    throw;

  Poss *APoss = new Poss(node1, false);
  RealPSet *ASet = new RealPSet(APoss);
  node->m_poss->AddPSet(ASet,true,true);

  alignModes.clear();
  alignModesSrc.clear();


  DistType BType;
  MatchDistsAndFillInWithStar(cont->m_BIndices,
			      CType, cont->m_CIndices,
			      BType, alignModes, alignModesSrc);

  RedistNode *node2 = new RedistNode(BType, cont->GetInputNameStr(2), alignModes, alignModesSrc);
  node2->AddInput(node->Input(1),node->InputConnNum(1));

  if (BType == node2->InputDataType(0).GetEffectiveDist())
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



string DistContToLocalContStatASumScatter::GetType() const
{
  return "DistContToLocalContStatASumScatter";
}

bool DistContToLocalContStatASumScatter::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Contraction::GetClass())
    throw;
  const Contraction *cont = (Contraction*)node;
  if (!cont->InputDataType(0).GetDist().HasNoReped())
    return false;
  if (cont->m_contIndices.empty())
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

  const DataTypeInfo &aInfo = ((DLANode*)(AConn->m_n))->DataType(AConn->m_num);
  if (aInfo.HasPerm())
    throw;
  const DistType &AType = aInfo.GetDist();

  EntryList sumDims;
  DimSet sumSet;
  string::iterator iter = cont->m_contIndices.begin();
  for(; iter != cont->m_contIndices.end(); ++iter) {
    size_t loc = cont->m_AIndices.find(*iter);
    if (loc != string::npos) {
      DistEntry entry = AType.m_dists[loc];
      sumDims.push_back(entry);
      DimVec vec = entry.DistEntryDims();
      sumSet.insert(vec.begin(), vec.end());
    }
  }

  DimVec alignModes, alignModesSrc;

  DistType BType;
  MatchDistsAndFillInWithStar(cont->m_BIndices,
			      AType, cont->m_AIndices,
			      BType, alignModes, alignModesSrc);

  RedistNode *node1 = NULL;
  RealPSet *BSet = NULL;

  if (BType != node->InputDataType(1).GetEffectiveDist()) {
    node1 = new RedistNode(BType, cont->GetInputNameStr(0), alignModes, alignModesSrc);
    node1->AddInput(node->Input(1),node->InputConnNum(1));

    Poss *BPoss = new Poss(node1, false);
    BSet = new RealPSet(BPoss);
    node->m_poss->AddPSet(BSet,true,true);
  }

  
  DistType CType;

  const DataTypeInfo &cInfo = cont->InputDataType(2);
  if (cInfo.HasPerm())
    throw;
  const DistType &CDestType = cInfo.GetDist();

  alignModes.clear();
  alignModesSrc.clear();
  MatchDistsAndFillInWithStar(cont->m_CIndices,
			      AType, cont->m_AIndices, 
			      CType, alignModes, alignModesSrc);

  TempVarNode *temp = new TempVarNode(CType, sumDims);
  temp->AddInput(node->Input(2),node->InputConnNum(2));
  
  Contraction *LCont = new Contraction(m_toLayer,  cont->m_alpha, COEFVALZERO, cont->m_type, 
				       cont->m_AIndices, cont->m_BIndices, cont->m_CIndices+cont->m_contIndices, cont->m_contIndices);
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
  for(int i = 0; iter != cont->m_CIndices.end(); ++iter, ++i) {
    size_t loc = cont->m_AIndices.find(*iter);
    if (loc != string::npos) {
      CType.m_dists[i] = AType.m_dists[loc]; 
      if (CType.m_dists[i] != CDestType.m_dists[i])
        notFinalType = true;
      DimVec dest = CDestType.m_dists[i].DistEntryDims();
      DimVecIter destIter = dest.begin();
      for(; destIter != dest.end(); ++destIter) {
	Dim dim = *destIter;
	DimSetIter find = sumSet.find(dim);
	if (find != sumSet.end()) {
	  CType.m_dists[i].AppendDim(dim);
	  sumSet.erase(find);
	}
      }
    }
    else {
      DimVec CDims;
      DimVec dest = CDestType.m_dists[i].DistEntryDims();
      DimVecIter destIter = dest.begin();
      for(; destIter != dest.end(); ++destIter) {
	Dim dim = *destIter;
	DimSetIter find = sumSet.find(dim);
	if (find != sumSet.end()) {
	  CDims.push_back(dim);
	  sumSet.erase(find);
	}
      }
      CType.m_dists[i].DimsToDistEntry(CDims);
    }
  }

//  if (!sumSet.empty())
//    throw;


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


    DimVec ident;
    IdentDimVec(CDestType.m_numDims, ident);
    
    RedistNode *finalRedist = new RedistNode(CDestType, node->GetInputNameStr(2), ident, ident);
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
  if (!cont->InputDataType(1).GetDist().HasNoReped())
    return false;
  if (cont->m_contIndices.empty())
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

  const DataTypeInfo &bInfo = ((DLANode*)(BConn->m_n))->DataType(BConn->m_num);
  if (bInfo.HasPerm())
    throw;
  const DistType &BType = bInfo.GetDist();

  EntryList sumDims;
  DimSet sumSet;
  string::iterator iter = cont->m_contIndices.begin();
  for(; iter != cont->m_contIndices.end(); ++iter) {
    size_t loc = cont->m_BIndices.find(*iter);
    if (loc != string::npos) {
      DistEntry entry = BType.m_dists[loc];
      sumDims.push_back(entry);
      DimVec vec = entry.DistEntryDims();
      sumSet.insert(vec.begin(), vec.end());
    }
  }

  DimVec alignModes, alignModesSrc;

  DistType AType;
  MatchDistsAndFillInWithStar(cont->m_AIndices,
			      BType, cont->m_BIndices,
			      AType, alignModes, alignModesSrc);

  RedistNode *node1 = NULL;
  RealPSet *ASet = NULL;

  const DataTypeInfo &AInfo = node->InputDataType(0);
  if (AInfo.HasPerm())
    throw;
  if (AType != AInfo.GetDist()) {
    node1 = new RedistNode(AType, cont->GetInputNameStr(1), alignModes, alignModesSrc);
    node1->AddInput(node->Input(0),node->InputConnNum(0));

    Poss *APoss = new Poss(node1, false);
    ASet = new RealPSet(APoss);
    node->m_poss->AddPSet(ASet,true,true);
  }

  //  cout << "AType " << AType.PrettyStr() << endl;

  
  DistType CType;

  const DataTypeInfo &cInfo = cont->InputDataType(2);
  if (cInfo.HasPerm())
    throw;
  const DistType &CDestType = cInfo.GetDist();

  alignModes.clear();
  alignModesSrc.clear();
  MatchDistsAndFillInWithStar(cont->m_CIndices,
			      BType, cont->m_BIndices, 
			      CType, alignModes, alignModesSrc);

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
  for(int i = 0; iter != cont->m_CIndices.end(); ++iter, ++i) {
    size_t loc = cont->m_BIndices.find(*iter);
    if (loc != string::npos) {
      CType.m_dists[i] = BType.m_dists[loc]; 
      if (CType.m_dists[i] != CDestType.m_dists[i])
        notFinalType = true;
      DimVec dest = CDestType.m_dists[i].DistEntryDims();
      DimVecIter destIter = dest.begin();
      for(; destIter != dest.end(); ++destIter) {
	Dim dim = *destIter;
	DimSetIter find = sumSet.find(dim);
	if (find != sumSet.end()) {
	  CType.m_dists[i].AppendDim(dim);
	  sumSet.erase(find);
	}
      }
    }
    else {
      DimVec CDims;
      DimVec dest = CDestType.m_dists[i].DistEntryDims();
      DimVecIter destIter = dest.begin();
      for(; destIter != dest.end(); ++destIter) {
	Dim dim = *destIter;
	DimSetIter find = sumSet.find(dim);
	if (find != sumSet.end()) {
	  CDims.push_back(dim);
	  sumSet.erase(find);
	}
      }
      CType.m_dists[i].DimsToDistEntry(CDims);
    }
  }

//  if (!sumSet.empty())
//    throw;

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

    DimVec ident;
    IdentDimVec(CDestType.m_numDims, ident);
    
    RedistNode *finalRedist = new RedistNode(CDestType, cont->GetInputNameStr(2), ident, ident);
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
				 DistType &final, DimVec &alignModes, DimVec &alignModesSrc)
{
  final.PrepForNumDims(indices.length());
  for(unsigned int i = 0; i < indices.length(); ++i) {
    char index = indices[i];
    size_t loc = matchingIndices.find(index);
    if (loc == string::npos) {
      final.m_dists[i].SetToStar();
    }
    else {
      alignModes.push_back(i);
      alignModesSrc.push_back(loc);
      final.m_dists[i] = matchingDists.m_dists[loc];
    }
  }
}


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
    CTun->SetUpStats(FULLUP,FULLUP,
		     NOTUP,NOTUP);
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
    CTun->SetAllStats(PARTUP);
  }
  if (isCIndex)
    CTun->SetIndepIters();
  
  Contraction *newCont = new Contraction(m_toLayer, cont->m_alpha, 
					 isCIndex ? cont->m_beta : COEFONE, 
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



bool ContractionLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Contraction::GetClass()) {
    const Contraction *cont = (Contraction*)node;
    if (cont->GetLayer() != m_fromLayer)
      return false;
    for (Dim dim = 0; dim < cont->InputNumDims(0); ++dim) {
      if (!(*(cont->InputLen(0,dim)) <= m_bs))
	return false;
    }
    for (Dim dim = 0; dim < cont->InputNumDims(2); ++dim) {
      if (!(*(cont->InputLen(2,dim)) <= m_bs))
	return false;
    }
    return true;
  }
  return false;
}

void ContractionLowerLayer::Apply(Node *node) const
{
  Contraction *cont = (Contraction*)node;
  cont->SetLayer(m_toLayer);
}

string ContractionLowerLayer::GetType() const
{ 
  return "Contraction lower layer " + LayerNumToStr(m_fromLayer) 
  + " to " + LayerNumToStr(m_toLayer);
}


void UpdateWithPermutation(Contraction *cont, ConnNum contInput, Permutation &perm)
{
  Node *inNode = cont->Input(contInput);
  ConnNum inNum = cont->InputConnNum(contInput);

  if (inNode->GetNodeClass() == RedistNode::GetClass()) {
    RedistNode *oldRedist = (RedistNode*)inNode;
    if (oldRedist->m_info.HasPerm())
      throw;
    RedistNode *newInput = new RedistNode(oldRedist->m_info.GetDist(), perm, 
					  oldRedist->m_align, oldRedist->m_alignModes, oldRedist->m_alignModesSrc);
    cont->m_poss->AddNode(newInput);
    newInput->AddInput(inNode->Input(0),inNode->InputConnNum(0));
    cont->ChangeInput2Way(inNode, inNum, newInput, 0);
    if (oldRedist->m_children.empty()) {
      cont->m_poss->DeleteChildAndCleanUp(oldRedist);
    }
  }
  else if (inNode->GetClass() == Permute::GetClass()) {
    Permute *oldPerm = (Permute*)inNode;
    Permutation composedPerm = oldPerm->m_permutation.ComposeWith(perm);
    Permute *newInput = new Permute(composedPerm, oldPerm->GetLayer());
    newInput->AddInput(inNode->Input(0),inNode->InputConnNum(0));
    cont->m_poss->AddNode(newInput);
    cont->ChangeInput2Way(inNode, inNum, newInput, 0);
    if (oldPerm->m_children.empty()) {
      cont->m_poss->DeleteChildAndCleanUp(oldPerm);
    }
  }
  else {
    Permute *newInput = new Permute(perm, SMLAYER);
    newInput->AddInput(inNode,inNum);
    cont->m_poss->AddNode(newInput);
    cont->ChangeInput2Way(inNode, inNum, newInput, 0);
  }
}

bool PermuteWhileUnpacking::CanApply(const Node *node) const
{
  const Contraction *cont = (Contraction*)node;
  if (m_type > 2)
    throw;
  if (!cont->m_needsPacking)
    return false;
  if (cont->Input(m_type)->GetNodeClass() == RedistNode::GetClass() ||
      cont->Input(m_type)->GetNodeClass() == Permute::GetClass())
    return true;
  return false;
}

void PermuteWhileUnpacking::Apply(Node *node) const
{
  Contraction *cont = (Contraction*)node;
  string newA, newB, newC;
  string inner;
  if (m_type == 0) {
    string aOuter;
    string bOuter;
    StringIter aIter = cont->m_AIndices.begin();
    for(; aIter != cont->m_AIndices.end(); ++aIter) {
      if (cont->m_contIndices.find(*aIter) != string::npos) {
	inner += (*aIter);
      }
      else {
	aOuter += (*aIter);
      }
    }
    newA = aOuter;
    newA += inner;
	
	
    StringIter bIter = cont->m_BIndices.begin();
    for(; bIter != cont->m_BIndices.end(); ++bIter) {
      if (inner.find(*bIter) == string::npos) {
	bOuter += (*bIter);
      }
    }
    newB = inner;
    newB += bOuter;

    newC = aOuter + bOuter;
  }
  else if (m_type == 1) {
    string aOuter, bOuter;
    StringIter bIter = cont->m_BIndices.begin();
    for(; bIter != cont->m_BIndices.end(); ++bIter) {
      if (cont->m_contIndices.find(*bIter) != string::npos) {
	inner += (*bIter);
      }
      else {
	bOuter += (*bIter);
      }
    }
    newB = inner;
    newB += bOuter;

    StringIter aIter = cont->m_AIndices.begin();
    for(; aIter != cont->m_AIndices.end(); ++aIter) {
      if (inner.find(*aIter) == string::npos) {
	aOuter += (*aIter);
      }
    }
    newA = aOuter;
    newA += inner;
	
    newC = aOuter + bOuter;
  }
  else if (m_type == 2) {
    newC = cont->m_CIndices;
    string aOuter, bOuter;

    StringIter cIter = cont->m_CIndices.begin();
    for(; cIter != cont->m_CIndices.end(); ++cIter) {
      if (cont->m_AIndices.find(*cIter) != string::npos) {	
	aOuter += (*cIter);
      } 
      else if (cont->m_BIndices.find(*cIter) != string::npos) {	
	bOuter += (*cIter);
      }
      else
	throw;
    }

    StringIter aIter = cont->m_AIndices.begin();
    for(; aIter != cont->m_AIndices.end(); ++aIter) {
      if (aOuter.find(*aIter) == string::npos) {
	inner += (*aIter);
      }
    }
    newA = aOuter;
    newA += inner;
    
    newB = inner;
    newB += bOuter;

    newC = aOuter + bOuter;
  }
  else
    throw;


  

      
  if (newA != cont->m_AIndices) {
    Permutation perm(cont->m_AIndices, newA);
    UpdateWithPermutation(cont, 0, perm);
  }

  if (newB != cont->m_BIndices) {
    Permutation perm(cont->m_BIndices, newB);
    UpdateWithPermutation(cont, 1, perm);
  }



  if (newC != cont->m_CIndices) {
    if (newC.size() < cont->InputDataType(2).GetDist().m_numDims)
      newC += inner;
  }

  if (newC != cont->m_CIndices)
  {

    Permutation perm(cont->m_CIndices, newC);
    UpdateWithPermutation(cont, 2, perm);
  }

  if (newC != cont->m_CIndices) {
    NodeConnVecIter iter = cont->m_children.begin();
    for(; iter != cont->m_children.end(); ++iter) {
      Node *child = (*iter)->m_n;
      if (child->GetNodeClass() != SumScatterUpdateNode::GetClass()
	  && child->GetNodeClass() != RedistNode::GetClass()
	  && child->GetNodeClass() != AllReduceNode::GetClass()) 
	{
	  Permute *newPermute = new Permute(newC, cont->m_CIndices, SMLAYER);
	  cont->m_poss->AddNode(newPermute);
	  cont->RedirectChildren(newPermute, 0);
	  newPermute->AddInput(cont, 0);
	}
    }
  }

  if (m_type == 0 || m_type == 1) {
    if (newC.size() < cont->m_CIndices.size())
      newC += cont->m_contIndices;

  }

  cont->m_AIndices = newA;
  cont->m_BIndices = newB;
  cont->m_CIndices = newC;
  cont->m_contIndices = inner;
  cont->m_needsPacking = false;
      
  return;
}



#endif


