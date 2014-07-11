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


#include "layers.h"

#if DOLLDLA

#include "loopSupport.h"
#include "gemm.h"

FullyUnrollLoop::FullyUnrollLoop(int maxNumIters) 
  : m_numIters(maxNumIters) 
{
  if (m_numIters <= 0)
    throw;
}

bool FullyUnrollLoop::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != SplitSingleIter::GetClass())
    throw;

  if (!node->IsPossTunnel(SETTUNIN))
    return false;

  const SplitSingleIter *split = (SplitSingleIter*)node;

  if (!split->m_isControlTun)
    return false;

  const Loop *loop = split->GetMyLoop();

  PossMMapConstIter iter = loop->m_posses.begin();
  for(; iter != loop->m_posses.end(); ++iter) {
    const Poss *poss = iter->second;
    PSetVecConstIter iter2 = poss->m_sets.begin();
    for( ; iter2 != poss->m_sets.end(); ++iter2) {
      if ((*iter2)->IsLoop())
	return false;
    }
  }
  
  unsigned int numExecs = split->NumberOfLoopExecs();
  if (!numExecs)
    throw;
  for(unsigned int i = 0; i < numExecs; ++i) {
    unsigned int numIters = split->NumIters(i);
    if (numIters != m_numIters)
      return false;
  }
  return true;
}

void RedirectChildrenThatArentPossTuns(Node *input, 
					  unsigned int num,
					  Node *newInput,
					  unsigned int newNum)
{
  for(unsigned int i = 0; i < input->m_children.size(); ++i) {
    NodeConn *output = input->m_children[i];
    if (output->m_num == num &&
	!output->m_n->IsPossTunnel()) 
      {
	//	cout << "redirecting child " << output->m_n->GetType() 
	//	     << " of " << input->GetType() << " to "
	//	     << newInput->GetType() << endl;
	input->RedirectChild(i, newInput, newNum);
	--i;
      }
  }
}

Poss* UnrollPoss(Poss *poss, Loop *loop, int numIters)
{
  NodeMap** map = new NodeMap*[numIters];
  Poss** newPosses = new Poss*[numIters];
  for (int dupNum = 0; dupNum < numIters; ++dupNum) {
    newPosses[dupNum] = new Poss;
    map[dupNum] = new NodeMap;
    newPosses[dupNum]->Duplicate(poss, *(map[dupNum]), true);
    newPosses[dupNum]->PatchAfterDuplicate(*(map[dupNum]), true);
  }
  
  Poss *rootPoss = newPosses[0];

  unsigned int numIn = poss->m_inTuns.size();
  //  cout << "orig numIn " << numIn << endl;
  for(unsigned int tunNum = 0; tunNum < numIn; ++tunNum) {
    //    cout << "***Working on tunNum " << tunNum << endl;
    LoopTunnel *possTunIn = (LoopTunnel*)(poss->m_inTuns[tunNum]);
    PossTunnel *newTun = NULL;
    PossTunnel *newOutTun = NULL;
    PossTunnel *outTun = NULL;

    if (!possTunIn->IsSplit()) {
      newTun = new PossTunnel(POSSTUNIN);

      Node *prev = NULL;
      Node *inputToOutTun = NULL;
      unsigned int inputNumToOutTun = 99999;
      
      outTun = possTunIn->GetMatchingOutTun();
      if (outTun) {
	if (!outTun->IsPossTunnel(POSSTUNOUT))
	  throw;
	if (outTun->IsCombine())
	  throw;

	inputToOutTun = outTun->Input(0);
	inputNumToOutTun = outTun->InputConnNum(0);
	//	cout << "tunNum " << tunNum << endl;
	//	cout << "matching outTun on rootPoss " << (*(map[0]))[outTun] << endl;
	//	cout << "input is " << inputToOutTun->GetNodeClass() << endl;
	//	cout << "input is " << inputToOutTun->GetType() << endl;
	/*
	if (inputToOutTun == possTunIn)
	  cout << "is same\n";
	else
	  cout << "is not same\n";
	*/
      }

      if (!outTun || (inputToOutTun == possTunIn)) {
	for(int dupNum = 0; dupNum < numIters; ++dupNum) {
	  Poss *dup = newPosses[dupNum];
	  LoopTunnel *possTunIn = (LoopTunnel*)(dup->m_inTuns[tunNum]);
	  RedirectChildrenThatArentPossTuns(possTunIn, 0, newTun, 0);
	  possTunIn->RemoveAllChildren2Way();
	}
      }
      else {
	prev = NULL;
	for(int dupNum = 0; dupNum < numIters; ++dupNum) {
	  Poss *dup = newPosses[dupNum];
	  LoopTunnel *possTunIn = (LoopTunnel*)(dup->m_inTuns[tunNum]);
	  if (dupNum) {
	    RedirectChildrenThatArentPossTuns(possTunIn, 0, prev, inputNumToOutTun);
	  }
	  else {
	    RedirectChildrenThatArentPossTuns(possTunIn, 0, newTun, 0);
	  }
	  //	  RedirectChildrenThatArentPossTuns(possTunIn, 0, newTun, 0);
	  NodeMapIter find = map[dupNum]->find(inputToOutTun);
	  if (find == map[dupNum]->end())
	    throw;
	  else {
	    prev = find->second;
	  }
	}
      }
      
      if (outTun) {
	newOutTun = new PossTunnel (POSSTUNOUT);
	if (prev)
	  newOutTun->AddInput(prev, inputNumToOutTun);
	else
	  newOutTun->AddInput(newTun, 0);
      }
    }
    else { // tunnel is a split
      newTun = new PossTunnel(POSSTUNIN);

      if (possTunIn->GetNodeClass() != SplitSingleIter::GetClass())
	throw;
      SplitSingleIter *split = (SplitSingleIter*)possTunIn;

      if (split->m_dir != PARTDOWN &&
	  split->m_dir != PARTRIGHT)
	throw;
      
      ViewMultipleIters *view = new ViewMultipleIters(split->m_dir,
						      split->GetMyLoop()->m_bsSize,
						      numIters);

      CombineMultipleIters *com = new CombineMultipleIters(split->m_dir,
							   split->GetMyLoop()->m_bsSize,
							   numIters);

      view->AddInput(newTun, 0);
      rootPoss->AddNode(view);
      rootPoss->AddNode(com);

      Node *inputToOutTun = NULL;
      unsigned int inputNumToOutTun = 99999;

      outTun = split->GetMatchingOutTun();
      if (!outTun)
	throw;
      if (outTun->GetNodeClass() != CombineSingleIter::GetClass())
	throw;
      if (!outTun->IsPossTunnel(POSSTUNOUT))
	throw;

      inputToOutTun = outTun->Input(1);
      inputNumToOutTun = outTun->InputConnNum(1);

      bool passThrough = inputToOutTun->IsPossTunnel(POSSTUNIN);
      
      for(int dupNum = 0; dupNum < numIters; ++dupNum) {
	Poss *dup = newPosses[dupNum];
	LoopTunnel *possTunIn = (LoopTunnel*)(dup->m_inTuns[tunNum]);

	//	cout << "old " << split << endl;
	/*
	cout << "orig possTunIn " << possTunIn << endl;
	if (passThrough)
	  cout << "pasthrough\n";
	else {
	  cout << "!pasthrough\n";
	}
	*/


	RedirectChildrenThatArentPossTuns(possTunIn, 1, view, dupNum);



	NodeConnVecIter iter = possTunIn->m_children.begin();
	for( ; iter != possTunIn->m_children.end(); ++iter) {
	  //part down and part right have 3 outputs for the split
	  // the 3rd (0-based) is the entire input
	  //this checks that all children only used the 1st (0-based)
	  // output 
	  if ((*iter)->m_num != 3) {
	    if (!(*iter)->m_n->IsPossTunnel(POSSTUNOUT)) {
	      cout << "num " << (*iter)->m_num << endl;
	      throw;
	    }
	  }
	}
	//This means it's a split and combine of read-only data
	// so it's just passed through
	//That's handled below
	if (!passThrough) {
	  NodeMapIter find = map[dupNum]->find(inputToOutTun);
	  if (find == map[dupNum]->end())
	    throw;
	  else {
	    com->AddInput(find->second, inputNumToOutTun);
	  }
	}
      }
      
      if (passThrough) {
	for (int i = 0; i < numIters+1; ++i)
	  com->AddInput(view, i);
      }
      else {
	com->AddInput(view, numIters);
      }

      newOutTun = new PossTunnel (POSSTUNOUT);
      newOutTun->AddInput(com, 0);
    }

    rootPoss->AddNode(newTun);
    rootPoss->m_inTuns.push_back(newTun);

    if (newOutTun) {
      rootPoss->AddNode(newOutTun);
      NodeVecIter outIter = rootPoss->m_outTuns.begin();
      Node *mappedOutTun = (*(map[0]))[outTun];
      if (!mappedOutTun)
	throw;
      mappedOutTun->RemoveAllInputs2Way();
      for(; outIter != rootPoss->m_outTuns.end(); ++outIter) {
	if (*outIter == mappedOutTun)
	  break;
      }
      if (outIter == rootPoss->m_outTuns.end())
	throw;
      NodeVecIter inserted = rootPoss->m_outTuns.insert(outIter, newOutTun);
      inserted++;
      if (*inserted != mappedOutTun)
	throw;
      rootPoss->m_outTuns.erase(inserted);
      rootPoss->DeleteNode(mappedOutTun);
      //      newOutTun->AddInput(newTun, 1);
      for(int dupNum = 1; dupNum < numIters; ++dupNum) {
	NodeMapIter find = map[dupNum]->find(outTun);
	if (find == map[dupNum]->end())
	  throw;
	else {
	  Node *possOutTun = find->second;
	  possOutTun->RemoveAllInputs2Way();
	  newPosses[dupNum]->DeleteNode(possOutTun);
	}
      }
    }
    else if (outTun)
      throw;
  }

  for(int dupNum = 0; dupNum < numIters; ++dupNum) {
    Poss *dup = newPosses[dupNum];
    //Notice this loops over the original number of intunnels
    //there should now be twice this many in the rootPoss
    for(unsigned int tunNum = 0; tunNum < numIn; ++tunNum) {
      // 0 since the 0th is deleted each time
      Node *tun = dup->m_inTuns[0];
      if (!tun->m_children.empty()) {
	//sanity check
	NodeConnVecIter iter = tun->m_children.begin();
	for(; iter != tun->m_children.end(); ++iter) {
	  Node *child = (*iter)->m_n;
	  if (!child->IsPossTunnel(POSSTUNOUT)) {
	    cout << poss->m_sets.size() << endl;
	    cout << child->GetNodeClass() << endl;
	    cout << child->GetType() << endl;
	    cout << "using output " << (*iter)->m_num << endl;
	    cout << "child " << child << endl;
	    cout << "tun " << tun << endl;
	    throw;
	  }
	}
      }
      NodeConnVecIter iter = tun->m_inputs.begin();
      for(; iter != tun->m_inputs.end(); ++iter)
	delete *iter;
      tun->m_inputs.clear();
      dup->DeleteNode(tun);
    }
  }
  
  //Just a sanity check
  NodeVecIter outIter = rootPoss->m_outTuns.begin();
  for(; outIter != rootPoss->m_outTuns.end(); ++outIter) {
    PossTunnel *outTun = (PossTunnel*)(*outIter);
    if (outTun->IsLoopTunnel()) {
      cout << "outTun " << outTun << endl;
      cout << outTun->GetNodeClass() << endl;
      cout << outTun->GetType() << endl;
      throw;
    }
  }    

  delete map[0];

  for(int dupNum = 1; dupNum < numIters; ++dupNum) {
    delete map[dupNum];
    Poss *dup = newPosses[dupNum];
    NodeVecIter iter = dup->m_possNodes.begin();
    for(; iter != dup->m_possNodes.end(); ++iter) {
      Node *node = *iter;
      //Sanity check
      if (node->IsPossTunnel())
	throw;
      node->m_poss = NULL;
      rootPoss->AddNode(node);
    }
    dup->m_possNodes.clear();
    
    delete dup;
  }

  delete [] newPosses;
  delete [] map;

  return rootPoss;
}


void FullyUnrollLoop::Apply(Node *node) const
{
  if (node->GetNodeClass() != SplitSingleIter::GetClass())
    throw;
  
  SplitSingleIter *split = (SplitSingleIter*)node;
  
  Loop *loop = split->GetMyLoop();

  PossMMapIter possIter = loop->m_posses.begin();

  Poss *rootPoss = UnrollPoss(possIter->second, loop, m_numIters);
  ++possIter;

  for (unsigned int i = 0; i < rootPoss->m_inTuns.size(); ++i) {
    Node *possInTun = rootPoss->m_inTuns[i];
    Node *loopInTun = loop->m_inTuns[i];

    possInTun->AddInput(loopInTun->Input(0), loopInTun->InputConnNum(0));
  }

  PSet *newSet = new PSet(rootPoss);
  loop->m_ownerPoss->AddPSet(newSet);

  if (newSet->m_outTuns.size() != loop->m_outTuns.size()) {
    cout << newSet->m_outTuns.size() << endl;
    cout << loop->m_outTuns.size() << endl;
    cout << rootPoss->m_outTuns.size() << endl;
    throw;
  }

  for(unsigned int i = 0; i < loop->m_outTuns.size(); ++i) {
    PossTunnel *loopTun = (PossTunnel*)(loop->m_outTuns[i]);
    if (!loopTun->m_children.empty()) {
      PossTunnel *setTun = (PossTunnel*)(newSet->m_outTuns[i]);
      loopTun->RedirectChildren(setTun);
    }
  }

  for(; possIter != loop->m_posses.end(); ++possIter) {
    Poss *poss = possIter->second;
    Poss *newPoss = UnrollPoss(poss, loop, m_numIters);
    newSet->AddPoss(newPoss);
  }

  

  NodeVecIter inIter = loop->m_inTuns.begin();
  for( ; inIter != loop->m_inTuns.end(); ++inIter) {
    Node *in = *inIter;
    in->RemoveAllInputs2Way();
    loop->m_ownerPoss->DeleteNode(in);
  }
  loop->m_inTuns.clear();

  for(unsigned int i = 0; i < loop->m_outTuns.size(); ++i) {
    PossTunnel *loopTun = (PossTunnel*)(loop->m_outTuns[i]);
    if (!loopTun->m_children.empty())
      throw;
    loop->m_ownerPoss->DeleteNode(loopTun);
    
  }
  loop->m_outTuns.clear();

  loop->m_ownerPoss->RemoveFromSets(loop);

  delete loop;
}


ViewMultipleIters::~ViewMultipleIters()
{
  if (m_sizes) {
    delete m_sizes;
    m_sizes = NULL;
  }
}

void ViewMultipleIters::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (m_inputs.size() != 1)
      throw;
    Input(0)->Prop();
    if (m_numIters <= 1)
      throw;
  }
}

void ViewMultipleIters::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const ViewMultipleIters *view = (ViewMultipleIters*)orig;
  m_partDir = view->m_partDir;
  m_bs = view->m_bs;
  m_numIters = view->m_numIters;
}

void ViewMultipleIters::ClearDataTypeCache()
{
  if (m_sizes) {
    delete m_sizes;
    m_sizes = NULL;
  }
}

void ViewMultipleIters::BuildDataTypeCache()
{
  if (m_sizes)
    return;
  m_sizes = new Sizes;

  m_info = InputDataType(0);

  switch (m_partDir) {
  case (PARTDOWN):
    {
      m_sizes->AddRepeatedSizes(BSSizeToSize(m_bs), GetInputM(0)->NumSizes(), 1);
      m_info.m_numRowsVar = BSSizeToVarName(m_bs);
      break;
    }
  case (PARTRIGHT):
    {
      m_sizes->AddRepeatedSizes(BSSizeToSize(m_bs), GetInputN(0)->NumSizes(), 1);
      m_info.m_numColsVar = BSSizeToVarName(m_bs);
      break;
    }
  default:
    throw;
  }
}

const Sizes* ViewMultipleIters::GetM(unsigned int num) const
{
  if (num > m_numIters) 
    throw;
  if (num == m_numIters) {
    return GetInputM(0);
  }

  if (m_partDir == PARTDOWN) {
    return m_sizes;
  }
  else if (m_partDir == PARTRIGHT) {
    return GetInputM(0);
  }
  else
    throw;
}

const Sizes* ViewMultipleIters::GetN(unsigned int num) const
{
  if (num > m_numIters)
    throw;

  if (num == m_numIters) {
    return GetInputN(0);
  }

  if (m_partDir == PARTRIGHT) {
    return m_sizes;
  }
  else if (m_partDir == PARTDOWN) {
    return GetInputN(0);
  }
  else
    throw;
}

Name ViewMultipleIters::GetName(unsigned int num) const
{
  if (num > m_numIters)
    throw;

  Name name = GetInputName(0);

  if (num == m_numIters) {
    return name;
  }

  switch (m_partDir)
    {
    case (PARTDOWN):
      name.m_name += "_VertSplit_";
      break;
    case (PARTRIGHT):
      name.m_name += "_HorzSplit_";
      break;
    default:
      throw;
    }
  name.m_name += std::to_string(num);
  return name;
}

void ViewMultipleIters::PrintCode(IndStream &out)
{
  const DataTypeInfo &type = InputDataType(0);
  string inName = GetInputNameStr(0);
  for(int i = 0; i < m_numIters; ++i) {
    out.Indent();
    *out << inName;
    if (m_partDir == PARTDOWN)
      *out << "_VertSplit_";
    else if (m_partDir == PARTRIGHT)
      *out << "_HorzSplit_";
    else
      throw;
    *out << i << " = " << inName;

    if (i) {
      *out << " + ";
      if (m_partDir == PARTDOWN) {
	if (!IsUnitStride(type.m_rowStride))
	  *out << type.m_rowStrideVar << " * ";
	if (m_bs == USELLDLAMU)
	  *out << MU_VAR_NAME << ";\n";
	else if (m_bs == USELLDLA2MU)
	  *out << "2 * " << MU_VAR_NAME << ";\n";
	else if (m_bs == USELLDLA3MU)
	  *out << "3 * " << MU_VAR_NAME << ";\n";
	else
	  throw;
      }
      else if (m_partDir == PARTRIGHT) {
	if (!IsUnitStride(type.m_colStride))
	  *out << type.m_colStrideVar << " * ";
	if (m_bs == USELLDLAMU)
	  *out << MU_VAR_NAME << ";\n";
	else if (m_bs == USELLDLA2MU)
	  *out << "2 * " << MU_VAR_NAME << ";\n";
	else if (m_bs == USELLDLA3MU)
	  *out << "3 * " << MU_VAR_NAME << ";\n";
	else
	  throw;
      }
      else
	throw;
    }
    *out << ";\n";
  }
}

void ViewMultipleIters::AddVariables(VarSet &set) const
{
  string inName = GetInputNameStr(0);
  if (m_partDir == PARTDOWN)
    inName += "_VertSplit_";
  else if (m_partDir == PARTRIGHT)
    inName += "_HorzSplit_";
  else
    throw;
  for(int i = 0; i < m_numIters; ++i) {
    Var var(inName, i);
    set.insert(var);
  }
}

void ViewMultipleIters::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  WRITE(m_bs);
  WRITE(m_partDir);
  WRITE(m_numIters);
}

void ViewMultipleIters::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLANode::UnflattenCore(in, info);
  READ(m_bs);
  READ(m_partDir);
  READ(m_numIters);
}


const DataTypeInfo& ViewMultipleIters::DataType(unsigned int num) const
{
  return m_info;
}


void CombineMultipleIters::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (m_inputs.size() != m_numIters+1)
      throw;

    if (m_numIters <= 1)
      throw;

    for (int i = 0; i < m_numIters+1; ++i)
      Input(i)->Prop();

    Node *in = Input(m_numIters);
    if (in->GetNodeClass() != ViewMultipleIters::GetClass())
      throw;
    ViewMultipleIters *inView = (ViewMultipleIters*)in;
    if (inView->m_bs != m_bs 
	|| inView->m_numIters != m_numIters
	|| inView->m_partDir != m_partDir)
      throw;


    switch (m_partDir) {
    case(PARTDOWN):
      {
	Size bs = BSSizeToSize(m_bs);
	for (int i = 0; i < m_numIters; ++i) {
	  if (*GetInputM(i) != bs)
	    throw;
	}	  
	break;
      }
    case(PARTRIGHT):
      {
	Size bs = BSSizeToSize(m_bs);
	for (int i = 0; i < m_numIters; ++i) {
	  if (*GetInputN(i) != bs)
	    throw;
	}
	break;
      }
    default:
      throw;
    }
  }
}

void CombineMultipleIters::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const CombineMultipleIters *view = (CombineMultipleIters*)orig;
  m_partDir = view->m_partDir;
  m_bs = view->m_bs;
  m_numIters = view->m_numIters;
}


const Sizes* CombineMultipleIters::GetM(unsigned int num) const
{
  if (num > 0) 
    throw;
  return GetInputM(m_numIters);
}

const Sizes* CombineMultipleIters::GetN(unsigned int num) const
{
  if (num > 0) 
    throw;
  return GetInputN(m_numIters);
}

Name CombineMultipleIters::GetName(unsigned int num) const
{
  if (num > 0) 
    throw;
  return GetInputName(m_numIters);
}

void CombineMultipleIters::PrintCode(IndStream &out)
{

}

void CombineMultipleIters::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  WRITE(m_bs);
  WRITE(m_partDir);
  WRITE(m_numIters);
}

void CombineMultipleIters::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLANode::UnflattenCore(in, info);
  READ(m_bs);
  READ(m_partDir);
  READ(m_numIters);
}

const DataTypeInfo& CombineMultipleIters::DataType(unsigned int num) const
{
  return InputDataType(m_numIters);
}

#endif

