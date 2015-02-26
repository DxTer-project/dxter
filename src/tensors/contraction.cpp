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

void FillIn(const DistEntry &CDestEntry, DistEntry &CTypeEntry, DimSet &sumSet, EntryList &sumDims);

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
  Size aSize = ((DLANode*)(cont->Input(0)))->TotalNumberOfLocalElements(cont->InputConnNum(0));
  Size bSize = ((DLANode*)(cont->Input(1)))->TotalNumberOfLocalElements(cont->InputConnNum(1));
  return aSize + bSize;
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

  
  RedistNode *node1 = new RedistNode(AType, GetAlignmentSource(cont,2), alignModes, alignModesSrc);
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

  RedistNode *node2 = new RedistNode(BType, GetAlignmentSource(cont,2), alignModes, alignModesSrc);
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
  Size bSize = ((DLANode*)(cont->Input(1)))->TotalNumberOfLocalElements(cont->InputConnNum(1));
  Size cSize = ((DLANode*)(cont->Input(2)))->TotalNumberOfLocalElements(cont->InputConnNum(2));
  return bSize + cSize;
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
    node1 = new RedistNode(BType, GetAlignmentSource(cont,0), alignModes, alignModesSrc);
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
  temp->m_align = GetAlignmentSource(cont,0);
  temp->m_alignModes = alignModes;
  temp->m_alignModesSrc = alignModesSrc;
  
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
      FillIn(CDestType.m_dists[i], CType.m_dists[i], sumSet, sumDims);
      if (CType.m_dists[i] != CDestType.m_dists[i])
        notFinalType = true;
    }
    else {
      CType.m_dists[i].SetToStar();
      FillIn(CDestType.m_dists[i], CType.m_dists[i], sumSet, sumDims);
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


    RedistNode *finalRedist = new RedistNode(CDestType, GetAlignmentSource((DLANode*)node,2), ident, ident);
    finalRedist->AddInput(sumSet->OutTun(0),0);
    Poss *redistPoss = new Poss(finalRedist, false);
    RealPSet *redistSet = new RealPSet(redistPoss);
    node->m_poss->AddPSet(redistSet,true,true);

    temp2->m_align = GetAlignmentSource((DLANode*)node,2);
    temp2->m_alignModes = ident;
    temp2->m_alignModesSrc = ident;
    
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
  Size aSize = ((DLANode*)(cont->Input(0)))->TotalNumberOfLocalElements(cont->InputConnNum(0));
  Size cSize = ((DLANode*)(cont->Input(2)))->TotalNumberOfLocalElements(cont->InputConnNum(2));
  return aSize + cSize;
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
    node1 = new RedistNode(AType, GetAlignmentSource(cont,1), alignModes, alignModesSrc);
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
  temp->m_align = GetAlignmentSource(cont,1);
  temp->m_alignModes = alignModes;
  temp->m_alignModesSrc = alignModesSrc;
    

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
      FillIn(CDestType.m_dists[i], CType.m_dists[i], sumSet, sumDims);
      if (CType.m_dists[i] != CDestType.m_dists[i])
        notFinalType = true;
    }
    else {
      CType.m_dists[i].SetToStar();
      FillIn(CDestType.m_dists[i], CType.m_dists[i], sumSet, sumDims);
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

    RedistNode *finalRedist = new RedistNode(CDestType, GetAlignmentSource(cont,2), ident, ident);
    finalRedist->AddInput(sumSet->OutTun(0),0);
    Poss *redistPoss = new Poss(finalRedist, false);
    RealPSet *redistSet = new RealPSet(redistPoss);
    node->m_poss->AddPSet(redistSet,true,true);

    temp2->m_align = GetAlignmentSource(cont,2);
    temp2->m_alignModes = ident;
    temp2->m_alignModesSrc = ident;
    
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
  alignModes.clear();
  alignModesSrc.clear();
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
  CTun->SetAdditive();
  
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
    return cont->GetLayer() == m_fromLayer;
  }
  throw;
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


void UpdateWithPermutation(Node *node, ConnNum nodeInput, Permutation &perm)
{
  Node *inNode = node->Input(nodeInput);
  ConnNum inNum = node->InputConnNum(nodeInput);

  if (inNode->GetNodeClass() == RedistNode::GetClass()) {
    RedistNode *oldRedist = (RedistNode*)inNode;
    if (oldRedist->m_info.HasPerm())
      throw;
    RedistNode *newInput = new RedistNode(oldRedist->m_info.GetDist(), perm, 
					  oldRedist->m_align, oldRedist->m_alignModes, oldRedist->m_alignModesSrc);
    node->m_poss->AddNode(newInput);
    newInput->AddInput(inNode->Input(0),inNode->InputConnNum(0));
    node->ChangeInput2Way(inNode, inNum, newInput, 0);
    if (oldRedist->m_children.empty()) {
      node->m_poss->DeleteChildAndCleanUp(oldRedist);
    }
  }
  else if (inNode->GetClass() == Permute::GetClass()) {
    Permute *oldPerm = (Permute*)inNode;
    Permutation composedPerm = oldPerm->m_permutation.ComposeWith(perm);
    Permute *newInput = new Permute(composedPerm, oldPerm->GetLayer());
    newInput->AddInput(inNode->Input(0),inNode->InputConnNum(0));
    node->m_poss->AddNode(newInput);
    node->ChangeInput2Way(inNode, inNum, newInput, 0);
    if (oldPerm->m_children.empty()) {
      node->m_poss->DeleteChildAndCleanUp(oldPerm);
    }
  }
  else if (inNode->IsTunnel(SETTUNOUT) && !((Tunnel*)inNode)->m_pset->IsLoop()) {
    bool updatedOne = false;
    bool needExternal = false;
    Tunnel *tun = (Tunnel*)inNode;
    for(auto possTunOut : tun->m_inputs) {
      if (possTunOut->m_n->m_inputs.size() != 1) {
	throw;
      }
      Node *in = possTunOut->m_n->Input(0);
      if (in->GetNodeClass() == RedistNode::GetClass()) {
	UpdateWithPermutation(possTunOut->m_n, 0, perm);
	updatedOne = true;
      }
      else if (in->GetNodeClass() == Permute::GetClass()) {
	throw;
      }
      else {
	needExternal = true;
      }
      if (needExternal && updatedOne)
	throw;
      if (needExternal) {
	Permute *newInput = new Permute(perm, SMLAYER);
	newInput->AddInput(inNode, 0);
	inNode->m_poss->AddNode(newInput);
	node->ChangeInput2Way(inNode, inNum, newInput, 0);
      }
      else if (!updatedOne) 
	throw;
    }
  }
  else {
    Permute *newInput = new Permute(perm, SMLAYER);
    newInput->AddInput(inNode,inNum);
    node->m_poss->AddNode(newInput);
    node->ChangeInput2Way(inNode, inNum, newInput, 0);
  }
}


bool PermIsFree(Contraction *cont, ConnNum contInput)
{
  Node *inNode = cont->Input(contInput);
  //  ConnNum inNum = cont->InputConnNum(contInput);

  if (inNode->GetNodeClass() == RedistNode::GetClass()) {
    return inNode->m_children.size() == 1;
  }
  else if (inNode->IsTunnel(SETTUNOUT)) {
    Tunnel *tun = (Tunnel*)inNode;
    if (tun->m_children.size() != 1)
      return false;
    return tun->m_pset->HasRedist();
  }
  else if (inNode->GetClass() == Permute::GetClass()) {
    return inNode->m_children.size() == 1;
  }
  else {
    return false;
  }
}

bool PermuteWhileUnpacking::CanApply(const Node *node) const
{
  if (CurrPhase != PACKOPTPHASE)
    return false;
  const Contraction *cont = (Contraction*)node;
  if (!cont->m_needsPacking)
    return false;
  for(ConnNum num = 0; num < 3; ++num) {
    const Node *in = cont->Input(num);
    if (in->GetNodeClass() == RedistNode::GetClass() ||
	in->GetNodeClass() == Permute::GetClass() ||
	in->GetNodeClass() == TempVarNode::GetClass())
      return true;
    if (in->IsTunnel(SETTUNOUT)) {
      BasePSet *set = ((Tunnel*)in)->m_pset;
      if (set->HasPermutableOut(FindInTunVec(set->m_outTuns,(Tunnel*)in)))
	return true; 
    }
  }
  return false;
}



void PermuteWhileUnpacking::Apply(Node *node) const

{
  Contraction *cont = (Contraction*)node;

  if (!cont->m_needsPacking)
    throw;

  string newAForA, newBForA, newCForA;
  string newAForB, newBForB, newCForB;
  string newAForC, newBForC, newCForC;
  string innerForA, innerForB, innerForC;
  Size aSize = 0;
  Size bSize = 0;
  Size cSize = 0;
  {
    string aOuter;
    string bOuter;
    StringIter aIter = cont->m_AIndices.begin();
    for(; aIter != cont->m_AIndices.end(); ++aIter) {
      if (cont->m_contIndices.find(*aIter) != string::npos) {
	innerForA += (*aIter);
      }
      else {
	aOuter += (*aIter);
      }
    }
    newAForA = aOuter;
    newAForA += innerForA;
	
	
    StringIter bIter = cont->m_BIndices.begin();
    for(; bIter != cont->m_BIndices.end(); ++bIter) {
      if (innerForA.find(*bIter) == string::npos)
	{
	  bOuter += (*bIter);
	}
    }
    newBForA = innerForA;
    newBForA += bOuter;

    newCForA = aOuter + bOuter;
    if (newCForA != cont->m_CIndices) {
      if (newCForA.size() < cont->InputDataType(2).GetDist().m_numDims)
	newCForA += innerForA;
    }

    if (newAForA != cont->m_AIndices
	&& !PermIsFree(cont, 0)) 
      {
	aSize += cont->TotalNumberOfInputLocalElements(0);
      }
    if (newBForA != cont->m_BIndices
	&& !PermIsFree(cont, 1)) 
      {
	aSize += cont->TotalNumberOfInputLocalElements(1);
      }
    if (newCForA != cont->m_CIndices) {
      if (!PermIsFree(cont, 2)) 
	{
	  aSize += cont->TotalNumberOfInputLocalElements(2);
	}
      if (cont->m_children.size() != 1) {
	aSize += cont->TotalNumberOfInputLocalElements(2);
      }
      else {
	Node *child = cont->Child(0);
	if ((child->IsTunnel(SETTUNIN) && !((Tunnel*)child)->m_pset->HasPermutableIn(FindInTunVec(((Tunnel*)child)->m_pset->m_inTuns,(Tunnel*)child)))
	    || (child->GetNodeClass() != RedistNode::GetClass() && child->GetNodeClass() != Permute::GetClass())) 
	  {
	    aSize += cont->TotalNumberOfInputLocalElements(2);
	  }
      }
    }
  }
  {
    string aOuter, bOuter;
    StringIter bIter = cont->m_BIndices.begin();
    for(; bIter != cont->m_BIndices.end(); ++bIter) {
      if (cont->m_contIndices.find(*bIter) != string::npos) {
	innerForB += (*bIter);
      }
      else {
	bOuter += (*bIter);
      }
    }
    newBForB = innerForB;
    newBForB += bOuter;

    StringIter aIter = cont->m_AIndices.begin();
    for(; aIter != cont->m_AIndices.end(); ++aIter) {
      if (innerForB.find(*aIter) == string::npos) {
	aOuter += (*aIter);
      }
    }
    newAForB = aOuter;
    newAForB += innerForB;
	
    newCForB = aOuter + bOuter;
    if (newCForB != cont->m_CIndices) {
      if (newCForB.size() < cont->InputDataType(2).GetDist().m_numDims)
	newCForB += innerForB;
    }


    if (newAForB != cont->m_AIndices
	&& !PermIsFree(cont, 0)) 
      {
	bSize += cont->TotalNumberOfInputLocalElements(0);
      }
    if (newBForB != cont->m_BIndices
	&& !PermIsFree(cont, 1)) 
      {
	bSize += cont->TotalNumberOfInputLocalElements(1);
      }
    if (newCForB != cont->m_CIndices) {
      if (!PermIsFree(cont, 2)) 
	{
	  bSize += cont->TotalNumberOfInputLocalElements(2);
	}
      if (cont->m_children.size() != 1) {
	bSize += cont->TotalNumberOfInputLocalElements(2);
      }
      else {
	Node *child = cont->Child(0);
	if ((child->IsTunnel(SETTUNIN) && !((Tunnel*)child)->m_pset->HasPermutableIn(FindInTunVec(((Tunnel*)child)->m_pset->m_inTuns,(Tunnel*)child)))
	    || (child->GetNodeClass() != RedistNode::GetClass() && child->GetNodeClass() != Permute::GetClass())) 
	  {
	    bSize += cont->TotalNumberOfInputLocalElements(2);
	  }
      }
    }
  }
  {
    string aOuter, bOuter;

    StringIter aIter = cont->m_AIndices.begin();
    for(; aIter != cont->m_AIndices.end(); ++aIter) {
      if (cont->m_contIndices.find(*aIter) != string::npos) {
	innerForC += (*aIter);
      }
    }

    StringIter cIter = cont->m_CIndices.begin();
    for(; cIter != cont->m_CIndices.end(); ++cIter) {
      if (innerForC.find(*cIter) == string::npos) {
	if (cont->m_AIndices.find(*cIter) != string::npos) {	
	  aOuter += (*cIter);
	} 
	else if (cont->m_BIndices.find(*cIter) != string::npos) {	
	  bOuter += (*cIter);
	}
	else
	  throw;
      }
    }

    newAForC = aOuter;
    newAForC += innerForC;
    
    newBForC = innerForC;
    newBForC += bOuter;

    newCForC = aOuter + bOuter;
    if (newCForC != cont->m_CIndices) {
      if (newCForC.size() < cont->InputDataType(2).GetDist().m_numDims)
	newCForC += innerForC;
    }

    if (newAForC != cont->m_AIndices
	&& !PermIsFree(cont, 0)) 
      {
	cSize += cont->TotalNumberOfInputLocalElements(0);
      }
    if (newBForC != cont->m_BIndices
	&& !PermIsFree(cont, 1)) 
      {
	cSize += cont->TotalNumberOfInputLocalElements(1);
      }
    if (newCForC != cont->m_CIndices) {
      if (!PermIsFree(cont, 2)) 
	{
	  cSize += cont->TotalNumberOfInputLocalElements(2);
	}
      if (cont->m_children.size() != 1) {
	cSize += cont->TotalNumberOfInputLocalElements(2);
      }
      else {
	Node *child = cont->Child(0);
	if ((child->IsTunnel(SETTUNIN) && !((Tunnel*)child)->m_pset->HasPermutableIn(FindInTunVec(((Tunnel*)child)->m_pset->m_inTuns,(Tunnel*)child)))
	    || (child->GetNodeClass() != RedistNode::GetClass() && child->GetNodeClass() != Permute::GetClass())) 
	  {
	    cSize += cont->TotalNumberOfInputLocalElements(2);
	  }
      }
    }
  }

  string newA, newB, newC, inner;

  bool fixC = false;
  if (cSize <= aSize && cSize <= bSize) {
    newA = newAForC;
    newB = newBForC;
    newC = newCForC;
    inner = innerForC;
  }
  else if (aSize <= bSize) {
    newA = newAForA;
    newB = newBForA;
    newC = newCForA;
    inner = innerForA;
    fixC = true;
  }
  else {
    newA = newAForB;
    newB = newBForB;
    newC = newCForB;
    inner = innerForB;
    fixC = true;
  }
      
  if (newA != cont->m_AIndices) {
    Permutation perm(cont->m_AIndices, newA);
    UpdateWithPermutation(cont, 0, perm);
  }

  if (newB != cont->m_BIndices) {
    Permutation perm(cont->m_BIndices, newB);
    UpdateWithPermutation(cont, 1, perm);
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
	  && child->GetNodeClass() != AllReduceNode::GetClass()
	  && (!child->IsTunnel(SETTUNIN) || !((Tunnel*)child)->m_pset->HasPermutableIn(FindInTunVec(((Tunnel*)child)->m_pset->m_inTuns,(Tunnel*)child))))
	{
	  Permute *newPermute = new Permute(newC, cont->m_CIndices, SMLAYER);
	  cont->m_poss->AddNode(newPermute);
	  cont->RedirectChildren(newPermute, 0);
	  newPermute->AddInput(cont, 0);
	}
    }
  }

  if (fixC) {
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

//if the final C type contains some grid mode over
// which we are summing, add all of the modes that must be summed
// with it to the intermediate C type mode with that grid mode
void FillIn(const DistEntry &CDestEntry, DistEntry &CTypeEntry, DimSet &sumSet, EntryList &sumDims)
{
  bool didSomething = false;
  DimVec destTmp = CDestEntry.DistEntryDims();
  DimVec dest = CTypeEntry.DistEntryDims();
  DimVecIter destIter = destTmp.begin();
  for(; destIter != destTmp.end(); ++destIter) {
    Dim dim = *destIter;
    DimSetIter find = sumSet.find(dim);
    if (find != sumSet.end()) {
      EntryListIter iter = sumDims.begin();
      for (; iter != sumDims.end(); ++ iter) {
	DistEntry entry = *iter;
	if (entry.ContainsDim(dim)) {
	  DimVec sumDimEntry = entry.DistEntryDims();
	  dest.insert(dest.end(), sumDimEntry.begin(), sumDimEntry.end());
	  DimVecIter sumDimEntryIter = sumDimEntry.begin();
	  for(; sumDimEntryIter != sumDimEntry.end(); ++sumDimEntryIter)
	    sumSet.erase(*sumDimEntryIter);
	  didSomething = true;
	}
      }
    }
  }
  CTypeEntry.DimsToDistEntry(dest);
}


#endif


