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

#include "copy.h"
#include "pack.h"
#include "packToCopyAndZero.h"
#include "partition.h"
#include "recombine.h"
#include "setToZero.h"

#if DOLLDLA

PackToCopyAndZero::PackToCopyAndZero(Layer fromLayer, Layer toLayer) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

string PackToCopyAndZero::GetType() const {
  return "PackToCopyAndZero";
}

bool PackToCopyAndZero::CanApply(const Node* node) const {
  if (node->GetNodeClass() == Pack::GetClass()) {
    return true;
  }
  throw;
}

void PackToCopyAndZero::Apply(Node* node) const {
  Pack* pack = (Pack*) node;

  Partition* partition;
  if (pack->PackDir() == HORIZONTAL) {
    partition = new Partition(m_toLayer, pack->PackDir(), pack->PackN());
  } else {
    partition = new Partition(m_toLayer, pack->PackDir(), pack->PackM());
  }
  partition->AddInput(pack->Input(1), pack->InputConnNum(1));

  Copy* copy = new Copy(m_toLayer);
  copy->AddInputs(4,
		  pack->Input(0), pack->InputConnNum(0),
		  partition, 0);

  SetToZero* setZero = new SetToZero(m_toLayer);
  setZero->AddInput(partition, 1);
  
  Recombine* recombine = new Recombine(m_toLayer, pack->PackDir());
  recombine->AddInputs(6,
		       copy, 0,
		       setZero, 0,
		       pack->Input(0), pack->InputConnNum(0));

  pack->m_poss->AddNode(partition);
  pack->m_poss->AddNode(copy);
  pack->m_poss->AddNode(setZero);
  pack->m_poss->AddNode(recombine);
  pack->RedirectChildren(recombine, 0);
  pack->m_poss->DeleteChildAndCleanUp(pack);
}

#endif // DOLLDLA
