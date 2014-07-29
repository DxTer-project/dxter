/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2014, The University of Texas and Bryan Marker

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

GraphIter::GraphIter(Poss *poss)
{
  m_poss = NULL;
  m_hasPrinted = false;
  Init(poss);
}

GraphIter::GraphIter(const GraphIter &iter) 
{
  m_poss = NULL;
  m_hasPrinted = false;
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
  m_hasPrinted = false;
  m_setIters = new PossMMapIter[poss->m_sets.size()];
  m_subIters = new GraphIterPtr[poss->m_sets.size()];
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    m_setIters[i] = m_poss->m_sets[i]->GetPosses().begin();
    m_subIters[i] = new GraphIter(m_setIters[i]->second);
  }
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
  m_hasPrinted = rhs.m_hasPrinted;
  m_setIters = new PossMMapIter[m_poss->m_sets.size()];
  m_subIters = new GraphIterPtr[m_poss->m_sets.size()];
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    m_setIters[i] = rhs.m_setIters[i];
    m_subIters[i] = new GraphIter(*(rhs.m_subIters[i]));
  }
  return *this;
}

bool GraphIter::Increment()
{
  m_hasPrinted = false;
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

void GraphIter::AddCurrPossVars(VarSet &set) const
{
  NodeVecConstIter iter = m_poss->m_possNodes.begin();
  for(; iter != m_poss->m_possNodes.end(); ++iter) {
    (*iter)->AddVariables(set);
  }
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    m_subIters[i]->AddCurrPossVars(set);
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
    PossMMap map = m_poss->m_sets[i]->GetPosses();
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
  }
  
  return tot;

}

void GraphIter::PrintRoot(IndStream &out, GraphNum whichGraph, bool currOnly, BasePSet *owner)
{
  if (!currOnly)
    Init(m_poss);

  unsigned int numPSets = m_poss->m_sets.size();

  bool keepGoing = true;
  GraphNum graphNum = 1;

  while (keepGoing) {
    if (currOnly || whichGraph == 0 || whichGraph == graphNum) {
      if (!currOnly)
	*out << "/*** Algorithm " << graphNum << " ***" << endl;
      else
	*out << "/***\n";
      *out << "\tUnique Num: " << m_poss->m_num << endl;
      *out << "\tChild of: " << m_poss->m_parent << endl;
      *out << "\tResult of transformations:" << endl;
      TransVec transVec;
      GetCurrTransVec(transVec);
      TransVecConstIter transIter = transVec.begin();
      for( ; transIter != transVec.end(); ++transIter)
	*out << "\t" << (*transIter)->GetType() << endl;
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
      /*
	Need to fix the following that came from inside of Poss::Print
	instead of doing this, BLIS should add variables
	if (m_pset && m_poss->m_pset->IsLoop()
	  && ((Loop*)(m_pset))->GetType() == BLISLOOP)
	{
	  string loopLevel = out.LoopLevel(1);
	  string idx = "idx" + loopLevel;
	  string dimLen = "dimLen" + loopLevel;
	  string bs = "bs" + loopLevel;
	  out.Indent();
	  *out << "dim_t " << idx << ", " << dimLen << ", " << bs << ";\n";
	}
      */
      
      //This actualy sets some stuff so it can print
      if (!m_poss->CanPrint()) {
	cout << "couldn't print\n";
	throw;
      }
      
      NodeVecConstIter nodeIter = m_poss->m_inTuns.begin();
      for(; nodeIter != m_poss->m_inTuns.end(); ++nodeIter) {
	(*nodeIter)->Print(out, whichGraph, this);
      }
      bool hasPrinted = true;
      while(hasPrinted) {
	hasPrinted = false;
	NodeVecConstIter nodeIter = m_poss->m_possNodes.begin();
	for( ; nodeIter != m_poss->m_possNodes.end(); ++nodeIter) {
	  //Don't bring the poss out tunnels until the end
	  // so the repartitioning code all goes after the loop body
	  if (!(*nodeIter)->HasPrinted() && !(*nodeIter)->IsTunnel(POSSTUNOUT)) {
	    (*nodeIter)->Print(out, whichGraph, this);
	    hasPrinted |= (*nodeIter)->HasPrinted();
	  }
	}
	for(unsigned int i = 0; i < numPSets; ++i) {
	  if (!m_subIters[i]->m_hasPrinted &&
	      m_setIters[i]->second->CanPrint()) 
	    {
	      out.Indent();
	      *out << "//**** (out of " << m_poss->m_sets[i]->GetPosses().size() << ")\n";
	      m_poss->m_sets[i]->PrePrint(out,m_setIters[i]->second);
	      ++out;
	      RealPSet *real = m_poss->m_sets[i]->GetReal();
	      real->SetInTunsAsPrinted();
	      m_subIters[i]->Print(out, whichGraph, m_poss->m_sets[i]);
	      --out;
	      m_poss->m_sets[i]->PostPrint(out,m_setIters[i]->second);
	      out.Indent();
	      *out << "//****\n";
	      hasPrinted = true;
	    }
	}
      }
      
      nodeIter = m_poss->m_outTuns.begin();
      for(; nodeIter != m_poss->m_outTuns.end(); ++nodeIter) {
	(*nodeIter)->Print(out, whichGraph, this);
	(*nodeIter)->SetPrinted();
      }
      
      nodeIter = owner->m_outTuns.begin();
      for(; nodeIter != owner->m_outTuns.end(); ++nodeIter) {
	(*nodeIter)->Print(out, whichGraph, this);
	(*nodeIter)->SetPrinted();
      }
      *out << endl;
      
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
  m_hasPrinted = true;
}

void GraphIter::Print(IndStream &out, GraphNum &graphNum, BasePSet *owner)
{
  m_hasPrinted = true;
  m_poss->ClearPrintedFromGraph();
  unsigned int numPSets = m_poss->m_sets.size();
  
  NodeVecConstIter nodeIter = m_poss->m_inTuns.begin();
  for(; nodeIter != m_poss->m_inTuns.end(); ++nodeIter) {
    (*nodeIter)->Print(out, graphNum, this);

    if (!(*nodeIter)->HasPrinted()) {
      cout << "tunnel input " << (*nodeIter)->GetType()
	   << "hasn't printed even though he should have\n";
    }
  }
  out.Indent();
  *out << "//------------------------------------//\n" << endl;
  bool hasPrinted = true;
  while(hasPrinted) {
    hasPrinted = false;
    NodeVecConstIter nodeIter = m_poss->m_possNodes.begin();
    for( ; nodeIter != m_poss->m_possNodes.end(); ++nodeIter) {
      Node *node = *nodeIter;
      //Don't print the poss out tunnels until the end
      // so the repartitioning code all goes after the loop body
      if (!node->HasPrinted()
          && !node->IsTunnel(POSSTUNOUT)
          && !node->IsTunnel(SETTUNOUT)
          && node->CanPrintCode(this))
	{
	  (*nodeIter)->Print(out, graphNum, this);
	  hasPrinted |= (*nodeIter)->HasPrinted();
	}
    }
    for(unsigned int i = 0; i < numPSets; ++i) {
      if (!m_subIters[i]->m_hasPrinted && m_poss->m_sets[i]->CanPrint(this)) {
	BasePSet *set = m_poss->m_sets[i];
	out.Indent();
	*out << "//**** (out of " << m_poss->m_sets[i]->GetPosses().size() << ")\n";
	out.Indent();
	*out << "//**** ";
	if (set->IsReal())
	  *out << "Is real\n";
	else
	  *out << "Is a shadow\n";
	RealPSet *real = set->GetReal();
	//	if (real->IsLoop())
	//	  *out << "is loop\n";
	//	*out << "real " << real << " instead of " << m_poss->m_sets[i] << endl;
	real->PrePrint(out,m_setIters[i]->second);
	++out;
	real->SetInTunsAsPrinted();
	m_subIters[i]->Print(out, graphNum, set);
	--out;
	real->PostPrint(out,m_setIters[i]->second);
	out.Indent();
	*out << "//****\n";

        hasPrinted = true;
      }
    }
  }
  
  
  *out << endl;
  out.Indent();
  *out << "//------------------------------------//" << endl;
  
  nodeIter = m_poss->m_outTuns.begin();
  for(; nodeIter != m_poss->m_outTuns.end(); ++nodeIter) {
    (*nodeIter)->Print(out, graphNum, this);
    (*nodeIter)->SetPrinted();
  }
  
  nodeIter = owner->m_outTuns.begin();
  for(; nodeIter != owner->m_outTuns.end(); ++nodeIter) {
    (*nodeIter)->Print(out, graphNum, NULL);
    (*nodeIter)->SetPrinted();
  }
  *out << endl;
  
  bool bad = false;
  
  nodeIter = m_poss->m_inTuns.begin();
  for(; nodeIter != m_poss->m_inTuns.end(); ++nodeIter) {
    if (!(*nodeIter)->HasPrinted()) {
      cout << (*nodeIter)->GetType() << " hasn't printed\n";
      bad = true;
    }
  }
  
  for(unsigned int i = 0; i < numPSets; ++i) {
    if (!m_subIters[i]->m_hasPrinted) {
      cout << "set not printed\n";
      for(unsigned int j = 0; j < m_poss->m_sets[i]->m_inTuns.size(); ++j) {
        Node *tun = m_poss->m_sets[i]->m_inTuns[j];
        cout << "in tun " << tun << endl;
        if (!tun->CanPrintCode(this))
          cout << "can't print\n";
      }
      m_subIters[i]->m_poss->CanPrint();
      m_subIters[i]->m_poss->ForcePrint();
      bad = true;
    }
  }
  
  nodeIter = m_poss->m_possNodes.begin();
  for(; nodeIter != m_poss->m_possNodes.end(); ++nodeIter) {
    if (!(*nodeIter)->HasPrinted()) {
      cout << (*nodeIter)->GetType() << " " << *nodeIter << " hasn't printed\n";
      cout << "Inputs are\n";
      (*nodeIter)->PrintInputs();
      cout << "Is it possible that the node is read only but ReadOnly doesn't return true?\n\n";
      bad = true;
      //      PrintSetConnections();
    }
  }
  
  if (bad) {
    cout << this << " is bad\n";
    cout << "contains " << m_poss->m_sets.size() << " posses\n";
    throw;
  }
}
