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

#include "hoistLoadToRegs.h"

#if DOLLDLA

#include "contiguousLoad.h"
#include "regLoadStore.h"

bool HoistLoadToRegs::OnlyFoundLoad(const LoopTunnel* setTunIn) const {
  for(auto childConn : setTunIn->m_children) {
    const LoopTunnel *possTunIn = (LoopTunnel*)(childConn->m_n);
    if (PossChildrenAreOnlyLoads(possTunIn)) {
      return true;
    }
  }
  return false;
}

bool HoistLoadToRegs::PossChildrenAreOnlyLoads(const LoopTunnel* possTunIn) const {
  bool foundLoad = false;
  bool foundSomethingElse = false;
  for(auto tunChildConn : possTunIn->m_children) {
    if (!tunChildConn->m_n->IsTunnel(POSSTUNOUT)) {
      if (tunChildConn->m_n->GetNodeClass() == LoadToRegs::GetClass()) {
	//	cout << "Found load" << endl;
	foundLoad = true;
      }
      else {
	/*	cout << "Also found something else" << endl;
		cout << tunChildConn->m_n->GetNodeClass() << endl;*/
	foundSomethingElse = true;
      }
    }
  }
  return foundLoad && !foundSomethingElse;
}

bool HoistLoadToRegs::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != LoopTunnel::GetClass()) {
    throw;
  }

  auto setTunIn = static_cast<const LoopTunnel*>(node);

  if (setTunIn->m_tunType != SETTUNIN) {
    return false;
  }

  if (!setTunIn->m_pset->IsReal()) {
    throw;
  }

  return OnlyFoundLoad(setTunIn);
}

void HoistLoadToRegs::RemoveLoadsFromPosses(LoopTunnel* setTunIn) const {
  auto set = static_cast<RealPSet*>(setTunIn->m_pset);
  for(int possNum = 0; possNum < setTunIn->m_children.size(); ++possNum) {
    auto possTunIn = static_cast<LoopTunnel*>(setTunIn->Child(possNum));
    auto poss = possTunIn->m_poss;
    if (PossChildrenAreOnlyLoads(possTunIn)) {
      for(auto tunChildConn : possTunIn->m_children) {
	if (!tunChildConn->m_n->IsTunnel(POSSTUNOUT)) {
	  auto oldLoad = static_cast<LoadToRegs*>(tunChildConn->m_n);
	  oldLoad->RedirectChildren(possTunIn, 0);
	  poss->DeleteChildAndCleanUp(oldLoad);
	} else {
	  LoopTunnel *possTunOut = static_cast<LoopTunnel*>(tunChildConn->m_n);
	  for(auto inputConn : possTunOut->m_inputs) {
	    if (!inputConn->m_n->IsTunnel(POSSTUNIN)) {
	      throw;
	      //handle this case by removing any stores
	    }
	  }
	}
      }
    } else {
      set->RemoveAndDeletePoss(poss, true);
      --possNum;
    }
  }
}

void HoistLoadToRegs::AddOneOuterLoad(LoopTunnel* setTunIn) const {
  auto newLoad = new LoadToRegs();
  setTunIn->m_poss->AddNode(newLoad);
  newLoad->AddInput(setTunIn->Input(0), setTunIn->InputConnNum(0));
  setTunIn->ChangeInput2Way(setTunIn->Input(0), setTunIn->InputConnNum(0), newLoad, 0);
}

void HoistLoadToRegs::Apply(Node *node) const {
  cout << "Applying HoistLoadToRegs" << endl;
  auto setTunIn = static_cast<LoopTunnel*>(node);
  if (!setTunIn->m_pset->IsReal()) {
    throw;
  }

  RemoveLoadsFromPosses(setTunIn);
  AddOneOuterLoad(setTunIn);

  if (setTunIn->GetMatchingOutTun()->m_children.size()) {
    //handle this case by adding a store between the set tun out and its children
    throw;
  }    
}

#endif // DOLLDLA
