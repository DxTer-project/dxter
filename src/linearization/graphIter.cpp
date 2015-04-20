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



#include "graphIter.h"
#include <iomanip>
#include "transform.h"
#include "contraction.h"
#include "realLoop.h"
#include "splitSingleIter.h"

TransTreeNode::~TransTreeNode()
{
  for(auto child : m_children)
    delete child;
  m_children.clear();
}


GraphIter::GraphIter(Poss *poss)
{
  m_poss = NULL;
  Init(poss);
}

GraphIter::GraphIter(const GraphIter &iter) 
{
  m_poss = NULL;
  *this = iter;
}

GraphIter::~GraphIter()
{

  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    delete m_subIters[i];
  }
  delete [] m_setIters;
  delete [] m_subIters;
  m_poss = NULL;

}

void GraphIter::Init(Poss *poss) 
{
  if (m_poss) {
    for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
      delete m_subIters[i];
    }
    delete [] m_setIters;
    delete [] m_subIters;
  }  
  m_poss = poss;
  m_setIters = new PossMMapIter[poss->m_sets.size()];
  m_subIters = new GraphIterPtr[poss->m_sets.size()];
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    m_setIters[i] = m_poss->m_sets[i]->GetPosses().begin();
    m_subIters[i] = new GraphIter(m_setIters[i]->second);
  }
  m_cost = -1;
}

GraphIter& GraphIter::operator=(const GraphIter &rhs)
{
  if (m_poss) {
    for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
      delete m_subIters[i];
    }
    delete [] m_setIters;
    delete [] m_subIters;
  }
  m_poss = rhs.m_poss;
  m_setIters = new PossMMapIter[m_poss->m_sets.size()];
  m_subIters = new GraphIterPtr[m_poss->m_sets.size()];
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    m_setIters[i] = rhs.m_setIters[i];
    m_subIters[i] = new GraphIter(*(rhs.m_subIters[i]));
  }
  m_cost = rhs.m_cost;
  return *this;
}

bool GraphIter::Increment()
{
  m_cost = -1;
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    bool ret = m_subIters[i]->Increment();
    if (!ret) 
      return false;
    else {
      BasePSet *set = m_poss->m_sets[i];
      delete m_subIters[i];
      ++(m_setIters[i]);
      PossMMap &map = set->GetPosses();
      if (m_setIters[i] == map.end()) {
	m_setIters[i] = map.begin();
	m_subIters[i] = new GraphIter(m_setIters[i]->second);
      }
      else {
	m_subIters[i] = new GraphIter(m_setIters[i]->second);
	return false;
      }
    }
  }
  return true;
}


void GraphIter::GetCurrTransVec(TransVec &transVec)
{
  transVec.insert(transVec.end(),
		  m_poss->m_transVec.begin(),m_poss->m_transVec.end());
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    m_subIters[i]->GetCurrTransVec(transVec);
  }
}

void GraphIter::PrintCurrTransformationTree(IndStream &out, unsigned int level, string set)
{
  ++out;
  for (auto trans : m_poss->m_transVec) {
    out.Indent();
    *out << "set " << set << ", level " << level << ", " << trans->GetType() << endl;
  }
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    string currSet = set + "." + std::to_string(i);
    m_subIters[i]->PrintCurrTransformationTree(out, level+1, currSet);
  }
  --out;
}

void GraphIter::AddCurrPossVars(VarSet &set) const
{
  NodeVecConstIter iter = m_poss->m_possNodes.begin();
  for(; iter != m_poss->m_possNodes.end(); ++iter) {
    (*iter)->AddVariables(set);
  }
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
#if !DOLOOPS
    m_subIters[i]->AddCurrPossVars(set);
#else
    BasePSet *base = m_poss->m_sets[i];
    RealPSet *real = base->GetReal();
    if (!real->IsLoop() ||
	!((RealLoop*)real)->IsUnrolled()) {
      m_subIters[i]->AddCurrPossVars(set);
    }
    else {
      RealLoop *loop = (RealLoop*)real;
      SplitSingleIter *con = (SplitSingleIter*)(loop->GetControl());
      int numIters = con->NumIters(0);
      for(int j = 0; j < numIters; ++j) {
	loop->SetCurrIter(j);
	m_subIters[i]->AddCurrPossVars(set);
      }
    }
#endif
  }
}


void GraphIter::EvalRoot(IndStream &out, GraphNum &graphNum, GraphNum whichGraph, GraphNum &optGraphs, double &optCosts)
{
  bool keepGoing = true;
  
  while (keepGoing) {
    if (whichGraph <= 0 || whichGraph == graphNum) {
      TransConstVec transList;
      Cost tot = Eval(transList);
#ifdef MATLAB
      *out << "cost(" << graphNum << ") = "
	   << setprecision(15) << tot << ";\n";
#else
      *out << "cost[" << graphNum << ",1] = "
	   << setprecision(15) << tot << ";\n";
#endif
      
#ifdef MATLAB
      *out << "refs(" << graphNum << ") = {[ ";
      TransConstVecIter iter = transList.begin();
      for(; iter != transList.end(); ++iter) {
        *out << (size_t)(*iter) << " ";
      }
      *out << "]};\n";
      if (!(graphNum % 1000)) {
        *out << "'loaded " << graphNum << "'\n";
      }
#endif
      if (optCosts <= 0 || tot < optCosts) {
        optCosts = tot;
        optGraphs = graphNum;
      }
    }
    
    ++graphNum;
    keepGoing = !Increment();
  }
}

//Update Eval()
Cost GraphIter::Eval(TransConstVec &transList)
{
  unsigned int numPSets = m_poss->m_sets.size();

  Cost tot = 0;
  
  TransVecIter iter2 = m_poss->m_transVec.begin();
  for(; iter2 != m_poss->m_transVec.end(); ++iter2) {
    transList.push_back(*iter2);
  }

  NodeVecConstIter iter = m_poss->m_possNodes.begin();
  for( ; iter != m_poss->m_possNodes.end(); ++iter) {
    DLANode *node =  (DLANode*)(*iter);
    double cost = node->GetCost();
    tot += cost;
  }

  for(unsigned int i = 0; i < numPSets; ++i) {
    tot += m_subIters[i]->Eval(transList);
  }

  m_cost = tot;
  
  return tot;
}

//Update Eval(TransConstVec &transList)
Cost GraphIter::Eval()
{
  unsigned int numPSets = m_poss->m_sets.size();

  if (m_cost != -1)
    return m_cost;

  Cost tot = 0;
  
  NodeVecConstIter iter = m_poss->m_possNodes.begin();
  for( ; iter != m_poss->m_possNodes.end(); ++iter) {
    DLANode *node =  (DLANode*)(*iter);
    double cost = node->GetCost();
    tot += cost;
  }

  for(unsigned int i = 0; i < numPSets; ++i) {
    tot += m_subIters[i]->Eval();
  }

  m_cost = tot;
  
  return tot;
}

Cost GraphIter::EvalAndSetBest()
{
  unsigned int numPSets = m_poss->m_sets.size();
  Cost tot = 0;
  
  NodeVecConstIter iter = m_poss->m_possNodes.begin();
  for( ; iter != m_poss->m_possNodes.end(); ++iter) {
    DLANode *node =  (DLANode*)(*iter);
    double cost = node->GetCost();
    tot += cost;
  }

  for(unsigned int i = 0; i < numPSets; ++i) {
    Cost optCost;
    PossMMap &map = m_poss->m_sets[i]->GetPosses();
    PossMMapIter iter = map.begin();
    m_subIters[i]->Init(iter->second);
    optCost = m_subIters[i]->EvalAndSetBest();
    ++iter;
    for( ; iter !=  map.end(); ++iter) {
      GraphIter tmp(iter->second);
      Cost tmpCost = tmp.EvalAndSetBest();
      if (tmpCost < optCost)  {
	*(m_subIters[i]) = tmp;
	m_setIters[i] = iter;
	optCost = tmpCost;
      }
    }
    tot += optCost;
  }

  m_cost = tot;
  
  return tot;

}

void GraphIter::PrintRoot(IndStream &out, GraphNum whichGraph, bool currOnly, BasePSet *owner)
{
  if (!currOnly)
    Init(m_poss);

  bool keepGoing = true;
  GraphNum graphNum = 1;

  while (keepGoing) {
    if (currOnly || whichGraph == 0 || whichGraph == graphNum) {
      if (!currOnly)
	*out << "/*** Algorithm " << graphNum << " ***" << endl;
      else if (whichGraph != 0)
	*out << "/*** Algorithm " << whichGraph << " ***" << endl;
      else
	*out << "/***\n";
      *out << "\tUnique Num: " << m_poss->m_num << endl;
      *out << "\tChild of: " << m_poss->m_parent << endl;
      *out << "\tResult of transformations:" << endl;
      PrintCurrTransformationTree(out, 0, "0");
      *out << "\n\tListed:\n";
      TransVec transVec;
      GetCurrTransVec(transVec);
      TransVecConstIter transIter = transVec.begin();
      for( ; transIter != transVec.end(); ++transIter)
	*out << "\t" << (*transIter)->GetType() << endl;
      *out << "\t\tCost = " << Eval() << endl;
      *out << "*****************************************/" << endl;
      
      VarSet set;
      AddCurrPossVars(set);
      VarSetIter varIter = set.begin();
#if DOTENSORS
      out.Indent();
      *out << "ObjShape tempShape;\n";
#endif
      for(; varIter != set.end(); ++varIter) {
	(*varIter).PrintDecl(out);
      }

      StrSet live;
      Print(out, NULL, live, 0);


      out.Indent();
      *out << "/*****************************************/" << endl;
      if (whichGraph != 0 || currOnly) {
	keepGoing = false;
      }
    }

    ++graphNum;
    if (keepGoing)
      keepGoing = !Increment();
  }
}

void GraphIter::Print(IndStream &out, BasePSet *owner, StrSet liveSet, Cost currCost)
{
  Linearizer lin(m_poss);
#if DOTENSORS
  lin.FindOptimalLinearization(liveSet);
  lin.InsertVecClearing(liveSet);
  //#if PRINTMEMCOSTS
  Cost highWater = -1;
  lin.m_lin.EnforceMemConstraint(currCost+lin.m_alwaysLiveCost, 1e10, liveSet, lin.m_alwaysLive, highWater);
  //#endif //PRINTMEMCOSTS

#else
  lin.FindAnyLinearization();
#endif

  for(auto in : m_poss->m_inTuns) {
    in->Print(out);
  }
  out.Indent();
  *out << "//------------------------------------//\n";
#if DOTENSORS
  out.Indent();
  *out << "\t////High water: " << highWater << endl;
#endif
  *out << endl;

  for(auto elem : lin.m_lin.m_order) {
#if PRINTMEMCOSTS
    cout << "currMemCost = " << elem->m_cost << endl;
#endif //PRINTMEMCOSTS
    if (!elem->IsSet()) {
      if (!elem->IsClear()) {
	*out << endl;
      }
      elem->Print(out);
    }
    else {
      SetLinElem *setElem = (SetLinElem*)elem;
      unsigned int num = FindInSetVec(m_poss->m_sets, setElem->m_set);
      GraphIter *subIter = m_subIters[num];
      Poss *poss = m_setIters[num]->second;
      
      setElem->Print(out, subIter, poss);
    }
  }
  
  
  *out << endl;
  out.Indent();
  *out << "//------------------------------------//" << endl;
  
  for(auto outTun : m_poss->m_outTuns) {
    outTun->Print(out);
  }
  if (owner) {
    for(auto outTun : owner->m_outTuns) {
      outTun->Print(out);
    }
  }

  *out << endl;
}




