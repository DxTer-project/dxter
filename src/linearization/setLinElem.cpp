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
#include "linearizer.h"
#include "realPSet.h"
#include "realLoop.h"
#include "splitSingleIter.h"

void SetLinElem::Print(IndStream &out)
{
  out.Indent();
  *out << "set for functionality "
       << m_set->GetFunctionalityString() << endl;
}

void SetLinElem::Print(IndStream &out, GraphIter *graphIter, Poss *poss)
{
  out.Indent();
  *out << "//**** (out of " << m_set->GetPosses().size() << ")\n";
  out.Indent();
  *out << "//**** ";
  if (m_set->IsReal()) {
    *out << "Is real\t" << ((RealPSet*)m_set)->m_shadows.size() << " shadows\n";
  }
  else {
    *out << "Is a shadow\t" << "of " << m_set->GetReal()->m_shadows.size() << " shadows\n";
  }

#if DOTENSORS
  out.Indent();
  *out << "\t//Outputs:\n";
  for(auto node : m_set->m_outTuns) {
    if (!node->m_children.empty()) {
      out.Indent();
      *out << "\t//  " << node->GetNameStr(0) << endl;
    }
  }
#endif //DOTENSORS


  RealPSet *real = m_set->GetReal();
  
#if !DOLOOPS
    real->PrePrint(out,poss);
    ++out;
    graphIter->Print(out, m_set, m_live, m_cost);
    --out;
    real->PostPrint(out,poss);
#else
  if (!real->IsLoop() ||
      !((RealLoop*)real)->IsUnrolled()) {
    real->PrePrint(out,poss);
    ++out;
    graphIter->Print(out, m_set, m_live, m_cost);
    --out;
    real->PostPrint(out,poss);
  }
  else {
    ++out;
    RealLoop *loop = (RealLoop*)real;
    SplitSingleIter *con = (SplitSingleIter*)(loop->GetControl());
    int numIters = con->NumIters(0);
    for(int j = 0; j < numIters; ++j) {
      loop->SetCurrIter(j);
      out.Indent();
      *out << "{\n";
      real->PrePrint(out,poss);
      graphIter->Print(out, m_set, m_live, m_cost);
      out.Indent();
      *out << "}\n";
    }
    --out;
  }
#endif

  out.Indent();
  *out << "//****\n";
}

//Used for determining cost, but not
// for when / where to clear variables - see Linearization::InsertVecClearing
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

StrSet SetLinElem::NewVars() const
{
  StrSet set;
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
	set.insert(outName);
      }
    }
  }
  return set;
}

VarCostMap SetLinElem::NewVarsAndCosts() const
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
#if DOTENSORS
	map[outName] = ((DLANode*)outTun)->MaxNumberOfLocalElements(0);
#elif TWOD
	map[outName] = ((DLANode*)outTun)->MaxNumberOfElements(0);
#else
	throw;
#endif
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

void SetLinElem::CacheLiveVars(const StrSet &stillLive)
{
  m_live = stillLive;
}

void SetLinElem::ClearCache()
{
  m_live.clear();
}
