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

#include "mvmulToSVMul.h"

#if DOLLDLA

#include "loopSupport.h"
#include "mvmul.h"
#include "svmul.h"
#include "vadd.h"

bool MVMulToSVMul::CanApply(const Node* node) const {
  if (node->GetNodeClass() == MVMul::GetClass()) {
    return true;
  }
  throw;
}

void MVMulToSVMul::Apply(Node* node) const {
  auto splitA = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitA->AddInput(node->Input(0), node->InputConnNum(0));
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();

  auto splitX = new SplitSingleIter(PARTDOWN, POSSTUNIN, false);
  splitX->AddInput(node->Input(1), node->InputConnNum(1));
  splitX->SetAllStats(FULLUP);
  splitX->SetIndepIters();

  auto tunY = new LoopTunnel(POSSTUNIN);
  tunY->AddInput(node->Input(2), node->InputConnNum(2));
  tunY->SetAllStats(PARTUP);

  auto svmul = new SVMul(m_toLayer);
  svmul->AddInputs(4,
		   splitX, 1,
		   splitA, 1);

  auto vadd = new VAdd(m_toLayer, COLVECTOR);
  vadd->AddInputs(4,
		  svmul, 0,
		  tunY, 0);

  auto comA = splitA->CreateMatchingCombine(0);
  auto comX = splitX->CreateMatchingCombine(0);

  auto outY = new LoopTunnel(POSSTUNOUT);
  outY->AddInput(vadd, 0);
  outY->AddInput(tunY, 1);
  outY->CopyTunnelInfo(tunY);

  Poss *loopPoss = new Poss(3, comA, comX, outY);
  RealLoop *loop = new RealLoop(LLDLALOOP, loopPoss, UnitBS);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(2), 0);
  node->m_poss->DeleteChildAndCleanUp(node);

  return;
}

#endif // DOLLDLA
