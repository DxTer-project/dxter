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


#include "linElem.h" 
#include "nodeLinElem.h"
#include "tempVarNode.h"

LinElem::LinElem() 
  : m_succ(NULL),
    m_addedToLinOrder(false)    
{
}

bool LinElem::FreeOfDataflowConstraints() const
{

  if (HasAdded())
    return false;
  for(auto input : m_inputs) {
    if (!input->HasAdded())
      return false;
  }
  for(auto pred : m_preds) {
    if (!pred->HasAdded())
      return false;
  }
  return true;
}

bool LinElem::CanAddToLinearOrder() const
{
  if (!FreeOfDataflowConstraints())
    return false;

  //if temp var node right before a set and there are other temp var nodes, 
  //  make sure they can all print before printing one since they should
  //  go right before the set

  if (IsNode()) {
    if (m_children.size() == 1) {
      const NodeLinElem *nodeElem = (NodeLinElem*)this;
      if (nodeElem->m_node->GetNodeClass() == TempVarNode::GetClass()) {
	if (m_succ)
	  throw;
	else {
	  const LinElem *child = m_children[0];

	  for (auto pred : child->m_preds) {
	    if (!pred->HasAdded()) {
	      return false;
	    }
	  }
	  for (auto input : child->m_inputs) {
	    if (!input->HasAdded()) {
	      if (input->IsNode()) {
		const NodeLinElem *nodeSetInput = (NodeLinElem*)input;
		if (nodeSetInput->m_node->GetNodeClass() == TempVarNode::GetClass()) {
		  if (!nodeSetInput->FreeOfDataflowConstraints())
		    return false;
		}
		else
		  return false;
	      }
	      else
		return false;
	    }
	  }
	}
      }
    }
  }
  return true;
}

void LinElem::AddInputIfUnique(LinElem *elem)
{
  for(auto input : m_inputs) {
    if (input == elem)
      return;
  }
  m_inputs.push_back(elem);
}

void LinElem::AddChildIfUnique(LinElem *elem)
{
  for(auto child : m_children) {
    if (child == elem)
      return;
  }
  m_children.push_back(elem);
}
