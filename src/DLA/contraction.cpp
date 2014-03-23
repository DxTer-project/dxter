#include "contraction.h"

#if DOTENSORS
#include "tensorRedist.h"

typedef vector<unsigned int *> ContType;
typedef ContType::iterator ContTypeIter;

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


Contraction::Contraction(Layer layer, Coef alpha, Coef beta, Type type, string indices)
:
  m_alpha(alpha),
  m_beta (beta),
  m_type (type),
  m_indices(indices)
{
  SetLayer(layer);
}

Node* Contraction::BlankInst()
{
  return new Contraction(ABSLAYER, COEFONE, COEFONE, REAL, "");
}

NodeType Contraction::GetType() const
{
  return "Contraction " + m_indices
    + LayerNumToStr(GetLayer());
}

void Contraction::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const Contraction *cont = (Contraction*)orig;
  m_alpha = cont->m_alpha;
  m_beta = cont->m_beta;
  m_type = cont->m_type;
  m_indices = cont->m_indices;
}

void Contraction::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_alpha);
  WRITE(m_beta);
  WRITE(m_type);
  out << m_indices << endl;
}

void Contraction::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  READ(m_alpha);
  READ(m_beta);
  READ(m_type);
  getline(in, m_indices);
}

const DistType& Contraction::GetDistType(unsigned int num) const
{
  return InputDistType(2);
}

void Contraction::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();
    m_cost = 0;
    cout << "improve Contraction::Prop code\n";
    cout << "reflect in DistContToLocalContStatC::RHSCostEstimate\n";
    DimVec dims = MapIndicesToDims(m_indices,GetInputName(0).m_indices);
    const Sizes *sizes = InputLocalLen(2,0);
    Dim numDims = InputNumDims(2);
    unsigned int totNumIters = sizes->NumSizes();
    for(unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
      Cost temp = 1;
      for (Dim dim = 0; dim < numDims; ++dim) {
	temp *= (*InputLocalLen(2,dim))[iteration];
      }
      DimVecConstIter iter = dims.begin();
      for(; iter != dims.end(); ++iter) {
	temp *= (*InputLocalLen(0,*iter))[iteration];
      }
      m_cost += temp;
    }
    m_cost *= 2;
  }
}

void Contraction::SanityCheck()
{
  
}

void Contraction::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "Contraction(";
  out << m_alpha;
  *out << ", " << GetInputName(0).str()
       << ", " << GetInputName(1).str()
       << ", ";
  out << m_beta;
  *out << ", " << GetInputName(2).str()
       << ", indices: " << m_indices
       << ");\n";
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

int DistContToLocalContStatC::CanApply(const Poss *poss, const Node *node, void **cache) const
{
  if (node->GetNodeClass() != Contraction::GetClass())
    throw;
  const Contraction *cont = (Contraction*)node;
  Dim numContDims = cont->m_indices.length();
  
  DimVec ADims = MapIndicesToDims(cont->m_indices,cont->GetInputName(0).m_indices);
  DimVec BDims = MapIndicesToDims(cont->m_indices,cont->GetInputName(1).m_indices);

  NodeConn *AConn = cont->InputConn(0);
  if (AConn->m_n->GetNodeClass() == RedistNode::GetClass())
    AConn = AConn->m_n->InputConn(0);
  NodeConn *BConn = cont->InputConn(1);
  if (BConn->m_n->GetNodeClass() == RedistNode::GetClass())
    BConn = BConn->m_n->InputConn(0);
  NodeConn *CConn = cont->InputConn(2);
  if (CConn->m_n->GetNodeClass() == RedistNode::GetClass())
    CConn = CConn->m_n->InputConn(0);



  const DistType &AType = ((DLANode*)(AConn->m_n))->GetDistType(AConn->m_num);
  const DistType &BType = ((DLANode*)(BConn->m_n))->GetDistType(BConn->m_num);
  const DistType &CType = ((DLANode*)(CConn->m_n))->GetDistType(CConn->m_num);

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
  /*
    dists ->  a numContDims-length c-style array of DimVecs;
              the kth DimVec holds the dimensions for
	      the k-th index's distribution
    thisDim -> this call is the thisDim^{th} index's call
    {A,B}Type -> DistType for the {A,B} tensor
    {A,B}Dims -> Map of contraction indices to dimensions of {A,B} (i.e.,
                 where those indices are in the tensor}
		 So the ADim[thisDim] distribution of AType is the
		 distribution of A for the current contraction index
    usedDims -> List of processes grid dimensions that have been 
                used for distribution up to this point in the recursion
    distOptions -> Vector of DistTypes 
  */
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

void DistContToLocalContStatC::Apply(Poss *poss, int num, Node *node, void **cache) const
{
  ContType *types = (ContType*)(*cache);
  unsigned int *entries = (*types)[num];

  Contraction *cont = (Contraction*)node;
  
  DimVec ADims = MapIndicesToDims(cont->m_indices,cont->GetInputName(0).m_indices);
  DimVec BDims = MapIndicesToDims(cont->m_indices,cont->GetInputName(1).m_indices);

  NodeConn *CConn = cont->InputConn(2);
  if (CConn->m_n->GetNodeClass() == RedistNode::GetClass())
    CConn = CConn->m_n->InputConn(0);

  const DistType &CType = ((DLANode*)(CConn->m_n))->GetDistType(CConn->m_num);

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
  RedistNode *node2 = new RedistNode(BType);
  RedistNode *node3 = new RedistNode(CType);
  Contraction *LCont = new Contraction(SMLAYER,  cont->m_alpha, cont->m_beta, cont->m_type, cont->m_indices);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node->Input(1),node->InputConnNum(1));
  node3->AddInput(node->Input(2),node->InputConnNum(2));
  LCont->AddInput(node1,0);
  LCont->AddInput(node2,0);
  LCont->AddInput(node3,0);

  bool sum = false;
  for (unsigned int i = 0; i < cont->m_indices.length() && !sum; ++i) {
    if (entries[i] != 0)
      sum = true;
  }

  DLANode *node4;
  if (sum)
    cout << "need different refinement code!\n";
  //  else {
    node4 = new RedistNode(cont->InputDistType(2));
    node4->AddInput(LCont, 0);
    //  }
    
  poss->AddNode(node1);
  poss->AddNode(node2);
  poss->AddNode(node3);
  poss->AddNode(LCont);
  poss->AddNode(node4);

  node->RedirectChildren(node4,0);
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


#endif
