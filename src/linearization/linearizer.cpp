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

#include "linearizer.h"
#include "poss.h"
#include "nodeLinElem.h"
#include "setLinElem.h"
#include "BasePSet.h"
#include "node.h"

Linearizer::Linearizer(const Poss *poss)
{
  PtrToLinElemMap map;
  for(auto node : poss->m_possNodes) {
    FindOrAdd(node, map);
  }
  for(auto set : poss->m_sets) {
    FindOrAdd(set, map);
  }
}


LinElem* Linearizer::FindOrAdd(const Node *node, PtrToLinElemMap &map)
{
  PtrToLinElemMapIter find = map.find(node);
  if (find != map.end())
    return find->second;

  if (node->IsTunnel()) {
    const Tunnel *tun = (Tunnel*)node;
    switch (tun->m_tunType)
      {
      case (POSSTUNIN):
      case (POSSTUNOUT):
	return NULL;
      case (SETTUNIN):
      case (SETTUNOUT):
	return FindOrAdd(tun->m_pset, map);
      default:
	throw;
      }
  }


  NodeLinElem *elem = new NodeLinElem(node);
  map[node] = elem;

  for(auto inputConn : node->m_inputs) {
    LinElem *inElem = FindOrAdd(inputConn->m_n, map);
    if (inElem)
      elem->AddInputIfUnique(inElem);
  }
  for(auto childConn : node->m_children) {
    LinElem *outElem = FindOrAdd(childConn->m_n, map);
    if (outElem)
      elem->AddChildIfUnique(outElem);
  }
  
  return elem;
}

LinElem* Linearizer::FindOrAdd(const BasePSet *set, PtrToLinElemMap &map)
{
  PtrToLinElemMapIter find = map.find(set);
  if (find != map.end())
    return find->second;
  
  SetLinElem *elem = new SetLinElem(set);
  map[set] = elem;

  for(auto inTun : set->m_inTuns) {
    for(auto inputConn : inTun->m_inputs) {
      LinElem *inElem = FindOrAdd(inputConn->m_n, map);
      if (inElem)
	elem->AddInputIfUnique(inElem);
    }
  }

  for(auto outTun : set->m_outTuns) {
    for(auto childConn : outTun->m_children) {
      LinElem *outElem = FindOrAdd(childConn->m_n, map);
      if (outElem)
	elem->AddChildIfUnique(outElem);
    }
  }
  
  return elem;
}
