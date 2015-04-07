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

#include "eliminateRecombinePartition.h"

#if DOLLDLA

#include "partition.h"
#include "recombine.h"

bool EliminateRecombinePartition::PartitionsAreIdentical(const Node* node) const {
  auto firstRecombine = static_cast<const Recombine*>(node);
  auto part = static_cast<const Partition*>(firstRecombine->Child(0));
  if (*part->GetM(0) == firstRecombine->GetInputNumRows(0)
      && *part->GetN(0) == firstRecombine->GetInputNumCols(0)
      && *part->GetM(1) == firstRecombine->GetInputNumRows(1)
      && *part->GetN(1) == firstRecombine->GetInputNumCols(1)) {
    cout << "Partitions are identical" << endl;
    return true;
  } else {
    //    cout << *part->GetM(0) << endl;
    cout << firstRecombine->GetInputNumRows(0) << endl;
    //    cout << *part->GetN(0) << endl;
    cout << firstRecombine->GetInputNumCols(0) << endl;
    //    cout << *part->GetM(1) << endl;
    cout << firstRecombine->GetInputNumRows(1) << endl;
    //    cout << *part->GetN(1) << endl;
    cout << firstRecombine->GetInputNumCols(1) << endl;
    return false;
  }
}

bool EliminateRecombinePartition::OutputIsSuperfluousPartition(const Node* node) const {
  if (node->NumChildrenOfOutput(0) == 2) {
    if (node->Child(0)->GetNodeClass() == Partition::GetClass()) {
      if (node->Child(1)->GetNodeClass() == Recombine::GetClass()) {
	auto er = node->Child(1);
	cout << "Child(1) name" << endl;
	cout << er->GetNodeClass() << endl;
	return PartitionsAreIdentical(node);
      } else {
	cout << "CANT APPLY" << endl;
	return false;
      }
    } else {
      cout << "CANT APPLY" << endl;
      return false;
    }
  } else {
    return false;
  }
}

bool EliminateRecombinePartition::CanApply(const Node* node) const {
  if (node->GetNodeClass() == Recombine::GetClass()) {
    return OutputIsSuperfluousPartition(node);
  }
  throw;
}

void EliminateRecombinePartition::Apply(Node* node) const {
  cout << "Applying EliminateRecombinePartition" << endl;
  auto superfluousPart = node->Child(0);
  auto er = node->Child(1);
  auto endRecombine = static_cast<Recombine*>(er);

  auto newEndRecombine = new Recombine(m_toLayer, endRecombine->GetDir());			     
  newEndRecombine->AddInputs(6,
			     endRecombine->Input(0), endRecombine->InputConnNum(0),
			     endRecombine->Input(1), endRecombine->InputConnNum(1),
			     node->Input(2), node->InputConnNum(2));

  endRecombine->RedirectChildren(newEndRecombine, 0);

  node->m_poss->DeleteChildAndCleanUp(endRecombine);

  superfluousPart->RedirectChildren(0, node->Input(0), node->InputConnNum(0));
  superfluousPart->RedirectChildren(1, node->Input(1), node->InputConnNum(1));

  node->m_poss->AddNode(newEndRecombine);
			    
  node->m_poss->DeleteChildAndCleanUp(superfluousPart);
  return;
}

#endif // DOLLDLA
