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

#include "copyColLoopRef.h"

#if DOLLDLA

#include "copy.h"
#include "loopSupport.h"

CopyColLoopRef::CopyColLoopRef(Layer fromLayer, Layer toLayer) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

bool CopyColLoopRef::CanApply(const Node* node) const {
  if (node->GetNodeClass() == Copy::GetClass()) {
    const Copy* copy = static_cast<const Copy*>(node);
    return *(copy->GetInputN(0)) > 1 &&
      *(copy->GetInputM(0)) > 1;
  }
  throw;
}

void CopyColLoopRef::Apply(Node* node) const {
  Copy* copy = static_cast<Copy*>(node);

  auto splitA = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitA->AddInput(copy->Input(0), copy->InputConnNum(0));
  splitA->SetUpStats(FULLUP, FULLUP,
		     FULLUP, FULLUP);
  splitA->SetIndepIters();

  auto splitB = new SplitSingleIter(PARTRIGHT, POSSTUNIN, false);
  splitB->AddInput(copy->Input(1), copy->InputConnNum(1));
  splitB->SetUpStats(FULLUP, NOTUP,
		     FULLUP, NOTUP);
  splitB->SetIndepIters();

  auto copyCol = new Copy(m_toLayer);
  copyCol->AddInput(splitA, 1);
  copyCol->AddInput(splitB, 1);

  auto comA = splitA->CreateMatchingCombine(0);
  auto comB = splitB->CreateMatchingCombine(1, 1, copyCol, 0);

  Poss* loopPoss = new Poss(2, comA, comB);
  RealLoop* loop = new RealLoop(LLDLALOOP, loopPoss, UnitBS);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

#endif // DOLLDLA
