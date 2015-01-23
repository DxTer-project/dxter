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


#include "setLinElem.h"
#include "tunnel.h"
#include "basePSet.h"

void SetLinElem::Print()
{
  throw;
}


StrVec SetLinElem::PossiblyDyingVars() const
{
  StrVec list;
  for(auto inTun : m_set->m_inTuns) {
    string inName = inTun->GetInputNameStr(0);
    bool alsoOutput = false;
    for(auto outTun : m_set->m_outTuns) {
      string outName = outTun->GetNameStr(0);
      if (inName == outName) {
	  if (outTun->m_children.size()) {
	    alsoOutput = true;
	    break;
	  }
      }
    }
    if (!alsoOutput)
      list.push_back(inName);
  }
  return list;
}

VarCostMap SetLinElem::NewVars() const
{
  VarCostMap map;
  for(auto outTun : m_set->m_outTuns) {
    if (outTun->m_children.size()) {
      string outName = outTun->GetNameStr(0);
      bool alsoInput = false;
      for(auto inTun : m_set->m_inTuns) {
	string inName = inTun->GetInputNameStr(0);
	if (inName == outName) {
	  alsoInput = true;
	  break;
	}
      }
      if (!alsoInput) {
	map[outName] = ((DLANode*)outTun)->MaxNumberOfLocalElements(0);
      }
    }
  }
  return map;
}

bool SetLinElem::UsesInputVar(const string &var) const
{
  for(auto inTun : m_set->m_inTuns) {
    if (inTun->GetInputNameStr(0) == var) {
      return true;
    }
  }
  return false;
}
