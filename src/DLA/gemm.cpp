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

#include "gemm.h"

#if DOBLIS||DOELEM||DOLLDLA
#include "elemRedist.h"
#include "string.h"
#include "helperNodes.h"
#include "loopSupport.h"
#include "blas.h"
#include "pack.h"

using namespace std;


#if DOELEM
bool IsDMGemm(const Node *node)
{
  if (node->GetNodeClass() != Gemm::GetClass())
    return false;
  const Gemm *gemm = (Gemm*)node;
  if (gemm->m_layer != DMLAYER)
    return false;
  else
    return true;
}
#endif



#if !DOLLDLA
Gemm::Gemm(Layer layer, Trans transA, Trans transB, Coef alpha, Coef beta, Type type)
: m_transA(transA),
m_transB(transB),
m_alpha(alpha),
m_beta(beta),
  m_type(type)
{
  SetLayer(layer);
#if DOBLIS
  m_comm = CORECOMM;
#endif
}
#else //DOLLDLA
Gemm::Gemm(Layer layer, Coef alpha, Coef beta, Type type)
: m_alpha(alpha),
m_beta(beta),
  m_type(type)
{
  SetLayer(layer);
#if DOBLIS
  m_comm = CORECOMM;
#endif
}
#endif //DOLLDLA


#if DOBLIS
bool Gemm::IsBLISParallelizable() const
{
  return GetLayer() == S3LAYER;
}

bool Gemm::IsParallel() const
{
  return m_comm != CORECOMM;
}

bool Gemm::RemoveParallelization()
{
  m_comm = CORECOMM;
  return false;
}

void Gemm::Parallelize(Comm comm)
{
  if (GetLayer() == S3LAYER)
    m_comm = comm;
  else
    throw;
}
#endif //DOBLIS


Node* Gemm::BlankInst()
{
  //These are default settings for a basic Gemm node
  return new Gemm(ABSLAYER, NORMAL, NORMAL, COEFONE, COEFONE, REAL);
}

NodeType Gemm::GetType() const
{
  return "Gemm "
  + TransToStr(m_transA)
  + TransToStr(m_transB) + " "
  + LayerNumToStr(GetLayer())
#if DOBLIS
    + " " + CommToStr(m_comm);
#else
  ;
#endif
}

void Gemm::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<3,1>::Duplicate(orig, shallow, possMerging);
  const Gemm *gemm = (Gemm*)orig;
#if !DOLLDLA
  m_transA = gemm->m_transA;
  m_transB = gemm->m_transB;
#endif
  m_alpha = gemm->m_alpha;
  m_beta = gemm->m_beta;
  m_type = gemm->m_type;
#if DOBLIS
  m_comm = gemm->m_comm;
#endif
}

void Gemm::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_transA);
  WRITE(m_transB);
  WRITE(m_alpha);
  WRITE(m_beta);
  WRITE(m_type);
#if DOBLIS
  WRITE(m_comm);
#endif
}

void Gemm::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  READ(m_transA);
  READ(m_transB);
  READ(m_alpha);
  READ(m_beta);
  READ(m_type);
#if DOBLIS
  READ(m_comm);
#endif
}

#if DOELEM
const DistType& Gemm::GetDistType(unsigned int num) const
{
#if DODPPHASE
  switch(GetLayer()) {
    case (ABSLAYER):
    case (DMLAYER):
      return MC_MR;
    case (SMLAYER):
      return InputDistType(2);
    default:
      throw;
  }
#else
  return InputDistType(2);
#endif
}
#endif

Phase Gemm::MaxPhase() const
{ 
  switch(GetLayer()) {
  case (ABSLAYER):
#if DODPPHASE
  case (DMLAYER):
    return DPPHASE;
  case (SMLAYER):
    return NUMPHASES;
#elif DOBLIS
    return SR1PHASE;
  case (S1LAYER):
    return SR2PHASE;
  case (S2LAYER):
    return SR3PHASE;
  case (S3LAYER):
    return NUMPHASES;
#endif
  default:
    throw;
  }
}

#if DOELEM
bool Gemm::DoNotCullDP() const
{
  return GetLayer() == DMLAYER;
}

bool Gemm::CanTransposeInputs() const
{
  return GetLayer() == SMLAYER;
}
#endif

Cost Gemm::GetCost(Layer layer, const Sizes *localDim1, const Sizes *localDim2, const Sizes *localDim3)
{
#if DOELEM
  if (layer == SMLAYER)
#elif DOBLIS
  if (layer == S1LAYER || layer == S2LAYER || layer == S3LAYER)
#else
  if (false)
#endif
    return TWO * GAMMA * localDim1->SumProds111(*localDim2, *localDim3);
  else
    throw;
}

void Gemm::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();
#if DOELEM
    if (GetLayer() == ABSLAYER || GetLayer() == DMLAYER) {
      if (InputDistType(2) != D_MC_MR) {
	cout << "input not D_MC_MR 7";
	throw;
      }
    }
#endif

    switch(GetLayer()) {
      case(ABSLAYER):
        m_cost = ZERO;
        break;
#if DOELEM
      case (DMLAYER):
        m_cost = ZERO;
        break;
      case (SMLAYER):
      {
        DLANode *in0 = (DLANode*)Input(0);
        unsigned int num0 = InputConnNum(0);
        DLANode *in2 = (DLANode*)Input(2);
        unsigned int num2 = InputConnNum(2);
        const Sizes *size1 = in0->LocalM(num0);
        const Sizes *size2 = in0->LocalN(num0);
        const Sizes *size3 = in2->LocalN(num2);
        m_cost = GetCost(SMLAYER, size1, size2, size3);
        break;
      }
#elif DOBLIS
      case (S1LAYER):
      case (S2LAYER):
      case (S3LAYER): 
	{
	  DLANode *in0 = (DLANode*)Input(0);
	  unsigned int num0 = InputConnNum(0);
	  DLANode *in2 = (DLANode*)Input(2);
	  unsigned int num2 = InputConnNum(2);
	  const Sizes *size1 = in0->LocalM(num0);
	  const Sizes *size2 = in0->LocalN(num0);
	  const Sizes *size3 = in2->LocalN(num2);
	  if (size1->NumSizes() != size2->NumSizes())
	    throw;
	  if (size2->NumSizes() != size3->NumSizes())
	    throw;
	  m_cost = GetCost(S3LAYER, size1, size2, size3) / NumCoresInComm(m_comm);
	  if (GetLayer() == S3LAYER && NumCoresInComm(m_comm) > 1) {
	    //BAM Should add cost for B panel, too.
	    Size numAElems = size1->SumProds11(*size2);
	    m_cost += AdditionalCostForBringingIntoL2(this, 0, numAElems, m_comm);
	  }
	  break;
	}
#endif
    default:
      throw;
    }
  }
}

#if DOELEM
void LocalGemmTransUpdate(DistType t0, DistType t1, Trans &transA, Trans &transB)
{
  if (transA == NORMAL) {
    if (t0 == D_MC_MR || t0 == D_MC_STAR || t0 == D_STAR_MC || t0 == D_STAR_MR)
      transA = NORMAL;
    else if (t0 == D_STAR_MC_H || t0 == D_STAR_MR_H || t0 == D_MR_STAR_H)
      transA = CONJTRANS;
    else if (t0 == D_STAR_MC_T || t0 == D_STAR_MR_T || t0 == D_MR_STAR_T)
      transA = TRANS;
    else if (t0 == D_STAR_STAR)
      transA = NORMAL;
    else {
      cout << "BAD 1!!!!!!\n";
      cout << TransToStr(transA) << endl;
      cout << DistTypeToStr(t0) << endl;
      throw;
    }
  }
  else if (transA == TRANS) {
    if (t0 == D_STAR_MC || t0 == D_STAR_MR || t0 == D_MC_MR || t0 == D_MC_STAR || t0 == D_MR_STAR)
      transA = TRANS;
    else if (t0 == D_MC_STAR_T || t0 == D_MR_STAR_T || t0 == D_STAR_MR_T || t0 == D_STAR_MC_T)
      transA = NORMAL;
    else if (t0 == D_STAR_STAR)
      transA = TRANS;
    else {
      cout << "BAD2!!!!!!\n";
      cout << DistTypeToStr(t0) << endl;
      throw;
    }
  }
  else if (transA == CONJTRANS) {
    if (t0 == D_STAR_MC || t0 == D_STAR_MR || t0 == D_MC_STAR || t0 == D_MC_MR)
      transA = CONJTRANS;
    else if (t0 == D_MC_STAR_H || t0 == D_MR_STAR_H || t0 == D_STAR_MC_H)
      transA = NORMAL;
    else if (t0 == D_STAR_STAR)
      transA = CONJTRANS;
    else {
      cout << "BAD 3!!!!!!\n";
      throw;
    }
  }
  else {
    cout << "BAD 4!!!!!!\n";
    throw;
  }
  if (t0 == D_STAR_STAR) {
    transB = transB;
  }
  else if (t0 == D_MC_MR) {
    if (transB == NORMAL) {
      if (t1 == D_MR_STAR || t1 == D_MC_STAR)
        transB = (NORMAL);
      else if (t1 == D_STAR_MR_T || t1 == D_STAR_MC_T)
        transB = (TRANS);
      else if (t1 == D_STAR_MR_H || t1 == D_STAR_MC_H)
        transB = (CONJTRANS);
      else {
        cout << "BAD 5!!!!!!\n";
        cout << DistTypeToStr(t1) << endl;
        throw;
      }
    }
    else if (transB == TRANS) {
      if (t1 == D_STAR_MR || t1 == D_STAR_MC)
        transB = (TRANS);
      else if (t1 == D_MR_STAR_T || t1 == D_MC_STAR_T)
        transB = (NORMAL);
      else {
        cout << "BAD 6!!!!!!\n";
        throw;
      }
    }
    else if (transB == CONJTRANS) {
      if (t1 == D_STAR_MR || t1 == D_STAR_MC)
        transB = (CONJTRANS);
      else if (t1 == D_MR_STAR_H || t1 == D_MC_STAR_H)
        transB = (NORMAL);
      else {
        cout << "BAD 7!!!!!!\n";
        throw;
      }
    }
    else  {
      cout << "BAD 8!!!!!!\n";
      throw;
    }
  }
  else if (transB == NORMAL) {
    if (t1 == D_STAR_MR || t1 == D_MC_MR)
      transB = (NORMAL);
    else if (t1 == D_MR_STAR_T)
      transB = (TRANS);
    else if (t1 == D_MR_STAR_H)
      transB = (CONJTRANS);
    else {
      cout << "BAD 9!!!!!!\n";
      throw;
    }
  }
  else if (transB == TRANS) {
    if (t1 == D_STAR_MR_T)
      transB = (NORMAL);
    else if (t1 == D_MR_STAR)
      transB = (TRANS);
    else if (t1 == D_MR_MC || t1 == D_MC_MR)
      transB = (TRANS);
    else {
      cout << "BAD 10!!!!!!\n";
      throw;
    }
  }
  else if (transB == CONJTRANS) {
    if (t1 == D_STAR_MR_H)
      transB = (NORMAL);
    else if (t1 == D_MR_STAR)
      transB = (CONJTRANS);
    else if (t1 == D_MR_MC || t1 == D_MC_MR)
      transB = (CONJTRANS);
    else {
      cout << "BAD 11!!!!!!" << DistTypeToStr(t0) << ", " << DistTypeToStr(t1) << "\n";
      throw;
    }
  }
  else
    cout << "BAD 12!!!!!!\n";
  
}
#endif

void Gemm::PrintCode(IndStream &out)
{
  out.Indent();
  if (GetLayer() == ABSLAYER ) {
    *out << "AbsGemm( "
	 << TransToStr(m_transA) << ", " << TransToStr(m_transB)
	 << ", \n\t";
    out << m_alpha;
    *out << ", "
    << GetInputName(0).str() << ", " << GetInputName(1).str()
    << ", ";
    out << m_beta;
    *out << ", " << GetInputName(2).str() << " );\n";
  }
#if DOELEM
  else if (GetLayer() == DMLAYER) {
    *out << "DistGemm( "
	 << TransToStr(m_transA) << ", " << TransToStr(m_transB)
    << ", \n\t";
    out << m_alpha;
    *out << ", "
    << GetInputName(0).str() << ", " << GetInputName(1).str()
    << ", ";
    out << m_beta;
    *out << ", " << GetInputName(2).str() << " );\n";
  }
  else if (GetLayer() == SMLAYER) {
    string transAStr, transBStr;
    DistType t0 = InputDistType(0);
    DistType t1 = InputDistType(1);
    Trans transA = m_transA;
    Trans transB = m_transB;
    LocalGemmTransUpdate(t0, t1, transA, transB);
    transAStr = TransToStr(transA);
    transBStr = TransToStr(transB);
    
    *out << "internal::LocalGemm( " << transAStr << ", " << transBStr << ", \n" << out.Tabs(1);
    out << m_alpha;
    *out << ","
    << GetInputName(0).str() << ", " << GetInputName(1).str() << ", \n" << out.Tabs(1);
    out << m_beta;
    *out << ", " << GetInputName(2).str() << " );\n";
  }
#elif DOBLIS
  else if (GetLayer() == S1LAYER ||
           GetLayer() == S2LAYER ||
           GetLayer() == S3LAYER) {
    string transAStr = TransToStr(m_transA);
    string transBStr = TransToStr(m_transB);
    
    if (GetLayer() == S1LAYER) {
      *out << "BLISGemmLimitedN( ";
      *out << transAStr << ", " << transBStr << ", \n" << out.Tabs(1);
      out << m_alpha;
      *out << ","
      << GetInputName(0).str() << ", " << GetInputName(1).str() << ", \n" << out.Tabs(1);
      out << m_beta;
      *out << ", " << GetInputName(2).str() << " );\n";
    }
    else if (GetLayer() == S2LAYER) {
      *out << "GemmRankKUpdate( ";
      *out << transAStr << ", " << transBStr << ", \n" << out.Tabs(1);
      out << m_alpha;
      *out << ","
      << GetInputName(0).str() << ", " << GetInputName(1).str() << ", \n" << out.Tabs(1);
      out << m_beta;
      *out << ", " << GetInputName(2).str() << " );\n";
    }
    else if (GetLayer() == S3LAYER) {
      if (m_comm == CORECOMM) 
	*out << "bli_gemm_ker_var2( ";
      else
	*out << "bli_gemm_ker_var2_par( ";
      out << m_alpha;
      *out<< ", &"
      << GetInputName(0).str() << ", &" << GetInputName(1).str() << ", \n"
      << out.Tabs(2);
      out << m_beta;
      *out << ", &" << GetInputName(2).str() << ", (gemm_t*)NULL";
      if (m_comm != CORECOMM)
	*out << ", " << CommToStr(GetSubComm(m_comm));
      *out << " );\n";
    }
    else
      throw;
    
  }
#endif
  else {
    throw;
  }
}

#if DOBLIS
void Gemm::UpdateInnerPackingMultiple(PackSize size)
{
  Node *input = Input(0);
  if (input->GetNodeClass() != Pack::GetClass())
    throw;
  Node *packInput = input->Input(1);
  if (packInput->GetNodeClass() != PackBuff::GetClass())
    throw;
  PackBuff *buff = (PackBuff*)packInput;
  if (buff->m_children.size() > 1)
    throw;
  buff->m_n = size;
}
#endif

#endif
