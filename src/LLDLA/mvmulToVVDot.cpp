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

#include "mvmulToVVDot.h"

#if DOLLDLA

#include "loopSupport.h"
#include "mvmul.h"
#include "vvdot.h"

bool MVMulToVVDot::CanApply(const Node* node) const {
  if (node->GetNodeClass() == MVMul::GetClass()) {
    return true;
  }
  throw;
}

void MVMulToVVDot::Apply(Node* node) const {
  auto splitA = new SplitSingleIter(PARTDOWN, POSSTUNIN, true);
  splitA->AddInput(node->Input(0), node->InputConnNum(0));
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();

  auto tunX = new LoopTunnel(POSSTUNIN);
  tunX->AddInput(node->Input(1), node->InputConnNum(1));
  tunX->SetAllStats(FULLUP);
  tunX->SetIndepIters();

  auto splitY = new SplitSingleIter(PARTDOWN, POSSTUNIN, false);
  splitY->AddInput(node->Input(2), node->InputConnNum(2));
  splitY->SetUpStats(FULLUP, FULLUP,
		     NOTUP, NOTUP);
  splitY->SetIndepIters();

  auto vvdot = new VVDot(m_toLayer);
  vvdot->AddInputs(6,
		   splitA, 0,
		   tunX, 0,
		   splitY, 0);

  auto comA = splitA->CreateMatchingCombine(0);
  auto comY = splitY->CreateMatchingCombine(1, 1, vvdot, 0);

  auto outX = new LoopTunnel(POSSTUNOUT);
  outX->AddInputs(4,
		  tunX, 0,
		  tunX, 1);
  outX->CopyTunnelInfo(tunX);

  Poss *loopPoss = new Poss(3, comA, outX, comY);
  RealLoop *loop = new RealLoop(LLDLALOOP, loopPoss, UnitBS);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(2), 0);
  node->m_poss->DeleteChildAndCleanUp(node);

  return;
}

#endif // DOLLDLA
