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



#include "lu.h"
#include "loopSupport.h"
#include "blas.h"

DistType LU::GetDistType(unsigned int num) const
{
  return InputDistType(num);
}
 
Phase LU::MaxPhase() const
{
  switch(GetLayer())
    {
    case (ABSLAYER):
#if DOSQM
      return SQR1PHASE;
#else
      throw;
#endif
      break;
    default:
      throw;
    }
}

void LU::SanityCheck()
{
  if (GetLayer() != ABSLAYER)
    throw;
  DLAOp<2,2>::SanityCheck();
}

void LU::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,2>::Prop();
    m_cost = ZERO;
  }
}

void LU::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "LU( " << GetInputNameStr(0) << ", " << GetInputNameStr(1) << " )\n";
}


string LULoopExp::GetType() const
{
  if (m_var == 5)
    return "LU Loop Exp var 5";
  else
    throw;
}

bool LULoopExp::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == LU::GetClass()) {
    const LU *lu = (LU*)node;
    if (lu->GetLayer() == m_fromLayer)
      return true;
  }
  return false;
}

void LULoopExp::Apply(Poss *poss, Node *node) const
{
  LU *lu = (LU*)node;
  Loop *loop;
  
  NodeConn *connA, *connB;
  connA = lu->m_inputs[0];
  connB = lu->m_inputs[1];
  
  switch(m_var) {
  case(5):
    loop = LUVar5Loop(connA->m_n, connA->m_num,
		      connB->m_n, connB->m_num,
		      m_toBLASLayer, m_toLAPACKLayer);
    break;
  default:
    throw;
  }
  
  poss->AddLoop(loop);
  
  node->RedirectChildren(0, loop->OutTun(0), 0);
  node->RedirectChildren(1, loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}


Loop* LUVar5Loop(Node *Ain, unsigned int Anum,
		 Node *Pin, unsigned int Pnum,
		 Layer BLASLayer, Layer LAPACKLayer)
{
  Split *splitA = new Split(PARTDIAG, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetUpStats(FULLUP, FULLUP,
		     PARTUP, PARTUP);
  
  Split *splitP = new Split(PARTDOWN, POSSTUNIN);
  splitP->AddInput(Pin, Pnum);
  splitP->SetUpStats(FULLUP, FULLUP,
                     NOTUP, NOTUP);

  Node *lu = new PanelLU(LAPACKLayer);
  lu->AddInputs(14, 
		splitA, 1,
		splitA, 2,
		splitA, 4,
		splitA, 5,
		splitA, 7,
		splitA, 8,
		splitP, 1);
		     

  Node *trxm = new Trxm(true, BLASLayer, LEFT, LOWER, UNIT, NORMAL, COEFONE, REAL);
  trxm->AddInputs(4,
		  lu, 2,
		  lu, 4);

  Node *gemm = new Gemm(BLASLayer, NORMAL, NORMAL, COEFNEGONE, COEFONE, REAL);
  gemm->AddInputs(6,
		  lu, 3,
		  trxm, 0,
		  lu, 5);
  
  Combine *comA = splitA->CreateMatchingCombine(6,
						1, lu, 0,
						2, lu, 1,
						4, lu, 2,
						5, lu, 3,
						7, trxm, 0,
						8, gemm, 0);
						
  
  Combine *comP = splitP->CreateMatchingCombine(1,
                                                1, lu, 6);
  
  Poss *loopPoss = new Poss(2, comA, comP);
  Loop *loop;
  if (BLASLayer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  else
    loop = new Loop(BLISLOOP, loopPoss, USEBLISMC);
  
  return loop;  
}


DistType PanelLU::GetDistType(unsigned int num) const
{
  return InputDistType(num);
}

Phase PanelLU::MaxPhase() const
{
  if (m_layer != S3LAYER)
    throw;
  return NUMPHASES;
}

void PanelLU::Prop()
{
  if (!IsValidCost(m_cost)) {
    m_cost = ZERO;
    DLAOp<7,7>::Prop();
  }
}

void PanelLU::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "PanelLU( &" 
       << "&" << GetInputNameStr(0) << ", &" << GetInputNameStr(1) << ",\n";
  out.Indent(1);
  *out << "&" << GetInputNameStr(2) << ", &" << GetInputNameStr(3) << ",\n";
  out.Indent(1);
  *out << "&" << GetInputNameStr(4) << ", &" << GetInputNameStr(5) << ",\n";
  out.Indent(1);
  *out << "&" << GetInputNameStr(6) << " );\n";
}
