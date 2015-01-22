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
      ClearPrintedRecursively();
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
      *out << "\t\tCost = " << m_poss->m_cost << endl;
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
#if PRINTEMPTY && (DOLLDLA == 0)
	    (*nodeIter)->PrintEmptyStatementIfOK(out);
#endif
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
	      BasePSet *set = m_poss->m_sets[i];
	      RealPSet *real = set->GetReal();
	      real->SetInTunsAsPrinted();
	      //Do this now so printing within here will properly empty variables
	      m_poss->m_sets[i]->m_flags |= SETHASPRINTEDFLAG;
	      m_subIters[i]->Print(out, whichGraph, m_poss->m_sets[i]);
	      --out;
	      m_poss->m_sets[i]->PostPrint(out,m_setIters[i]->second);

#if PRINTEMPTY && (DOLLDLA == 0)
	      NodeVecIter tunnelIter = m_poss->m_sets[i]->m_outTuns.begin();
	      for(; tunnelIter != m_poss->m_sets[i]->m_outTuns.end(); ++tunnelIter) {
		(*tunnelIter)->PrintEmptyStatementIfOK(out);
	      }
	      tunnelIter = m_poss->m_sets[i]->m_inTuns.begin();
	      for(; tunnelIter != m_poss->m_sets[i]->m_inTuns.end(); ++tunnelIter) {
		(*tunnelIter)->PrintEmptyStatementIfOK(out);
	      }
#endif

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

  for(auto in : m_poss->m_inTuns) {
    in->Print(out, graphNum, this);

    if (!in->HasPrinted()) {
      cout << "tunnel input " << in->GetType()
	   << "hasn't printed even though he should have\n";
    }
  }
  out.Indent();
  *out << "//------------------------------------//\n" << endl;
  bool hasPrinted = true;
  while(hasPrinted) {
    hasPrinted = false;
    for(auto node : m_poss->m_possNodes) {
      //Don't print the poss out tunnels until the end
      // so the repartitioning code all goes after the loop body
      if (!node->HasPrinted()
          && !node->IsTunnel(POSSTUNOUT)
          && !node->IsTunnel(SETTUNOUT)
          && node->CanPrintCode(this))
	{
	  node->Print(out, graphNum, this);
	  hasPrinted |= node->HasPrinted();
#if PRINTEMPTY && (DOLLDLA == 0)
	  node->PrintEmptyStatementIfOK(out);
#endif
	

	}
    }
    for(unsigned int i = 0; i < numPSets; ++i) {
      if (!m_subIters[i]->m_hasPrinted && m_poss->m_sets[i]->CanPrint(this)) {
	BasePSet *set = m_poss->m_sets[i];
	out.Indent();
	*out << "//**** (out of " << m_poss->m_sets[i]->GetPosses().size() << ")\n";
	out.Indent();
	*out << "//**** ";
	if (set->IsReal()) {
	  *out << "Is real\t" << ((RealPSet*)set)->m_shadows.size() << " shadows\n";
	}
	else {
	  *out << "Is a shadow\t" << "of " << set->GetReal()->m_shadows.size() << " shadows\n";
	}
#if DOTENSORS
	out.Indent();
	*out << "\t//Outputs:\n";
	for(auto node : set->m_outTuns) {
	  if (!node->m_children.empty()) {
	    out.Indent();
	    *out << "\t//  " << node->GetNameStr(0) << endl;
	  }
	}
#endif //DOTENSORS
	
	RealPSet *real = set->GetReal();
	set->m_flags |= SETHASPRINTEDFLAG;
	
	if (!real->IsLoop() ||
	    !((RealLoop*)real)->IsUnrolled()) {
	  real->PrePrint(out,m_setIters[i]->second);
	  ++out;
	  real->SetInTunsAsPrinted();
	  m_subIters[i]->Print(out, graphNum, set);
	  --out;
	  real->PostPrint(out,m_setIters[i]->second);
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
	    m_subIters[i]->ClearPrintedRecursively();
	    real->PrePrint(out,m_setIters[i]->second);
	    real->SetInTunsAsPrinted();
	    m_subIters[i]->Print(out, graphNum, set);
	    out.Indent();
	    *out << "}\n";
	  }
	  --out;
	}

#if PRINTEMPTY && (DOLLDLA == 0)
	for (auto outTun : m_poss->m_sets[i]->m_outTuns) {
	  outTun->PrintEmptyStatementIfOK(out);
	}
	for (auto inTun : m_poss->m_sets[i]->m_inTuns) { 
	  inTun->PrintEmptyStatementIfOK(out);
	}
#endif

	out.Indent();
	*out << "//****\n";

        hasPrinted = true;
      }
    }
  }
  
  
  *out << endl;
  out.Indent();
  *out << "//------------------------------------//" << endl;
  
  for(auto outTun : m_poss->m_outTuns) {
    outTun->Print(out, graphNum, this);
    outTun->SetPrinted();
  }
  
  for(auto outTun : owner->m_outTuns) {
    outTun->Print(out, graphNum, NULL);
    outTun->SetPrinted();
  }
  *out << endl;

  bool bad = false;
  
  for(auto inTun : m_poss->m_inTuns) {
    if (!inTun->HasPrinted()) {
      cout << inTun->GetType() << " hasn't printed\n";
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
  
  for(auto node : m_poss->m_possNodes) {
    if (!node->HasPrinted()) {
      cout << node->GetType() << " " << node << " hasn't printed\n";
      cout << "on " << node->m_poss << endl;
      cout << "Inputs are\n";
      node->PrintInputs();
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

void GraphIter::ClearPrintedRecursively()
{
  m_hasPrinted = false;
  m_poss->ClearPrintedFromGraph();
  unsigned int numPSets = m_poss->m_sets.size();
  for(unsigned int i = 0; i < numPSets; ++i) {
    m_poss->m_sets[i]->m_flags &= ~SETHASPRINTEDFLAG;
    m_subIters[i]->ClearPrintedRecursively();
  }
}


