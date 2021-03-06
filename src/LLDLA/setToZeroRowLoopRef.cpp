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

#include "loopSupport.h"
#include "setToZero.h"
#include "setToZeroRowLoopRef.h"

#if DOLLDLA

SetToZeroRowLoopRef::SetToZeroRowLoopRef(Layer fromLayer, Layer toLayer) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

bool SetToZeroRowLoopRef::CanApply(const Node* node) const {
  if (node->GetNodeClass() == SetToZero::GetClass()) {
    const SetToZero* setZero = static_cast<const SetToZero*>(node);
    return *(setZero->GetInputM(0)) > 1;
  }
  LOG_FAIL("replacement for throw call");
  throw;
}

void SetToZeroRowLoopRef::Apply(Node* node) const {
  SetToZero* setZero = static_cast<SetToZero*>(node);

  auto splitA = new SplitSingleIter(PARTDOWN, POSSTUNIN, true);
  splitA->AddInput(setZero->Input(0), setZero->InputConnNum(0));
  splitA->SetUpStats(FULLUP, FULLUP,
		     NOTUP, NOTUP);
  splitA->SetIndepIters();

  auto setRowZero = new SetToZero(m_toLayer);
  setRowZero->AddInput(splitA, 1);

  auto comA = splitA->CreateMatchingCombine(1, 1, setRowZero, 0);

  Poss* loopPoss = new Poss(1, comA);
  RealLoop* loop = new RealLoop(LLDLALOOP, loopPoss, UnitBS);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(0), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

#endif // DOLLDLA
