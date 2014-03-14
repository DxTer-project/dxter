#include "contraction.h"

#if DOTENSORS
#include "tensorRedist.h"

void RecursivelyFindDistributions(DimVec *dists, Dim thisDim, 
				  const DistType &AType, const DimVec &ADims,
				  const DistType &BType, const DimVec &BDims,
				  DimSet &usedDims,
				  DistTypeVec *distOptions);
void AddUnusedDimsForDistType(DimVec *dists, Dim numDims, 
			       DimSet &usedDims,
			      Dim numDimsUsed,
			       DistTypeVec *distOptions);

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
    unsigned int totNumIters = sizes->NumSizes();
    for(unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
      Cost temp = 1;
      Dim numDims = InputNumDims(2);
      for (Dim dim = 1; dim < numDims; ++dim) {
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
  throw;
}

void Contraction::PrintCode(IndStream &out)
{
  throw;
}

Phase Contraction::MaxPhase() const
{
  throw;
}

bool DistContToLocalContStatC::CanApply(const Poss *poss, const Node *node) const
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
  NodeConn *BConn = cont->InputConn(0);
   if (BConn->m_n->GetNodeClass() == RedistNode::GetClass())
     BConn = BConn->m_n->InputConn(0);
  NodeConn *CConn = cont->InputConn(0);
   if (CConn->m_n->GetNodeClass() == RedistNode::GetClass())
     CConn = CConn->m_n->InputConn(0);



  const DistType &AType = ((DLANode*)(AConn->m_n))->GetDistType(AConn->m_num);
  const DistType &BType = ((DLANode*)(BConn->m_n))->GetDistType(BConn->m_num);
  const DistType &CType = ((DLANode*)(CConn->m_n))->GetDistType(CConn->m_num);

  DimVec *dists = new DimVec[numContDims];
  
  DimSet usedDims = CType.UsedGridDims();

  DistTypeVec *distOptions = new DistTypeVec;

  RecursivelyFindDistributions(dists, 0, AType, ADims, BType, BDims, usedDims, distOptions);

  delete [] dists;
}

void RecursivelyFindDistributions(DimVec *dists, Dim thisDim, 
				  const DistType &AType, const DimVec &ADims,
				  const DistType &BType, const DimVec &BDims,
				  DimSet &usedDims,
				  DistTypeVec *distOptions)
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
    AddUnusedDimsForDistType(dists, ADims.size(), usedDims, usedDims.size(), distOptions);
    return;
  }
  //Fist, call recursively with * for this dim
  RecursivelyFindDistributions(dists, thisDim+1, 
			       AType, ADims,
			       BType, BDims,
			       usedDims,
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
				 usedDimsTemp,
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
				 usedDimsTemp,
				 distOptions);
  }
}

void AddUnusedDimsForDistType(DimVec *dists, Dim numDims, 
			       DimSet &usedDims,
			      Dim numDimsUsed,
			       DistTypeVec *distOptions)
{
  if (numDims == numDimsUsed) {
    distOptions->push_back(new DistType(dists));
    return;
  }
asdfasdflkjasdf
}

void DistContToLocalContStatC::Apply(Poss *poss, Node *node) const
{
  adsflkj
}

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


#endif
