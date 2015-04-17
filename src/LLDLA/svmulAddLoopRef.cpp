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

#include "svmulAddLoopRef.h"

#if DOLLDLA

#include "svmulAdd.h"

SVMulAddLoopRef::SVMulAddLoopRef(Layer fromLayer, Layer toLayer, VecType vtype)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_vtype = vtype;
}

string SVMulAddLoopRef::GetType() const
{
  switch(m_vtype)
    {
    case(ROWVECTOR):
      return "SVMulAddLoopRef - row vector";
    case(COLVECTOR):
      return "SVMulAddLoopRef - column vector";
    default:
      LOG_FAIL("replacement for throw call");
    }
}

bool SVMulAddLoopRef::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == SVMulAdd::GetClass()) {
    const SVMulAdd *svmulAdd = static_cast<const SVMulAdd*>(node);
    if (svmulAdd->GetLayer() != m_fromLayer) {
      return false;
    }
    if (m_vtype == ROWVECTOR) {
      return *(svmulAdd->GetInputN(1)) > svmulAdd->GetVecRegWidth() &&
	svmulAdd->GetInputN(1)->EvenlyDivisibleBy(svmulAdd->GetVecRegWidth());
    } 
    else if (m_vtype == COLVECTOR) {
      return *(svmulAdd->GetInputM(1)) > svmulAdd->GetVecRegWidth() &&
	svmulAdd->GetInputM(1)->EvenlyDivisibleBy(svmulAdd->GetVecRegWidth());
    } 
    else {
      LOG_FAIL("replacement for throw call");
      throw;
    }
  }

  cout << "ERROR: Cannot apply SVMulAddLoopRef to a non SVMulAdd node" << endl;
  LOG_FAIL("replacement for throw call");
  throw;
}

void SVMulAddLoopRef::Apply(Node *node) const {
  SVMulAdd* svmulAdd = static_cast<SVMulAdd*>(node);

  SplitSingleIter* splitA = new SplitSingleIter(m_vtype == COLVECTOR ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  splitA->AddInput(svmulAdd->Input(1), svmulAdd->InputConnNum(1));

  SplitSingleIter* splitB = new SplitSingleIter(m_vtype == COLVECTOR ? PARTDOWN : PARTRIGHT, POSSTUNIN, false);
  splitB->AddInput(svmulAdd->Input(2), svmulAdd->InputConnNum(2));

  if (m_vtype == COLVECTOR) {
    splitA->SetUpStats(FULLUP, FULLUP,
		       NOTUP, NOTUP);

    splitB->SetUpStats(FULLUP, FULLUP,
		       NOTUP, NOTUP);
  } else {
    splitA->SetUpStats(FULLUP, NOTUP,
		       FULLUP, NOTUP);

    splitB->SetUpStats(FULLUP, NOTUP,
		       FULLUP, NOTUP);
  }

  splitA->SetIndepIters();
  splitB->SetIndepIters();

  auto scalarTun = new LoopTunnel(POSSTUNIN);
  scalarTun->AddInput(svmulAdd->Input(0),svmulAdd->InputConnNum(0));
  scalarTun->SetAllStats(FULLUP);
  scalarTun->SetIndepIters();

  auto newMul = new SVMulAdd(svmulAdd->m_layer, svmulAdd->GetVecType());
  newMul->SetLayer(m_toLayer);
  newMul->AddInputs(6,
		    scalarTun, 0,
		    splitA, 1,
		    splitB, 1);

  auto scalarTunOut = new LoopTunnel(POSSTUNOUT);
  scalarTunOut->AddInput(scalarTun, 0);
  scalarTunOut->AddInput(scalarTun, 1);
  scalarTunOut->CopyTunnelInfo(scalarTun);

  auto comA = splitA->CreateMatchingCombine(0);

  auto comB = splitB->CreateMatchingCombine(1, 
					    1, newMul, 0);

  auto loopPoss = new Poss(3, scalarTunOut, comA, comB);

  RealLoop* loop;
  if (svmulAdd->GetDataType() == REAL_SINGLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuSingle);
  } else if (svmulAdd->GetDataType() == REAL_DOUBLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuDouble);
  } else {
    cout << "Error: Bad data type in vadd apply\n";
    LOG_FAIL("replacement for throw call");
    throw;
  }

  loop->SetDimName(m_vtype == COLVECTOR ? DIMM : DIMN);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(2), 0);
  node->m_poss->DeleteChildAndCleanUp(node);

}


#endif // DOLLDLA
