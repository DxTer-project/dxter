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


#include "layers.h"

#if DOLLDLA

#include "loopSupport.h"
#include "gemm.h"

FullyUnrollLoop::FullyUnrollLoop(int maxNumIters) 
  : m_numIters(maxNumIters) 
{
  if (m_numIters <= 0)
    LOG_FAIL("replacement for throw call");
}

bool FullyUnrollLoop::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != SplitSingleIter::GetClass()) {
    cout << "Error: Attempted to apply FullyUnrollLoop to non SplitSingleIter node\n";
    LOG_FAIL("replacement for throw call");
  }

  if (!node->IsTunnel(SETTUNIN))
    return false;

  const SplitSingleIter *split = (SplitSingleIter*)node;

  if (!split->m_isControlTun)
    return false;

  const BasePSet *loop = dynamic_cast<const BasePSet*>(split->GetMyLoop());
  
  unsigned int numExecs = split->NumberOfLoopExecs();
  if (!numExecs) {
    cout << "Error: Attempted to unroll loop with 0 executions\n";
    LOG_FAIL("replacement for throw call");
  }

  for(unsigned int i = 0; i < numExecs; ++i) {
    unsigned int numIters = split->NumIters(i);
    if (numIters != m_numIters)
      return false;
  }

  const PossMMap &map = loop->GetPosses();
  PossMMapConstIter iter = map.begin();
  for(; iter != map.end(); ++iter) {
    const Poss *poss = iter->second;
    if (poss->ContainsNonLoopCode())
      return true;
  }
  return false;
}

void RedirectChildrenThatArentTuns(Node *input, 
					  ConnNum num,
					  Node *newInput,
				       ConnNum newNum,
				       TunType type)
{
  for(unsigned int i = 0; i < input->m_children.size(); ++i) {
    NodeConn *output = input->m_children[i];
    if (output->m_num == num &&
	!output->m_n->IsTunnel(type)) 
      {
	//	cout << "redirecting child " << output->m_n->GetType() 
	//	     << " of " << input->GetType() << " to "
	//	     << newInput->GetType() << endl;
	input->RedirectChild(i, newInput, newNum);
	--i;
      }
  }
}

Poss* UnrollPoss(Poss *poss, LoopInterface *loop, int numIters)
{
  NodeMap** map = new NodeMap*[numIters];
  Poss** newPosses = new Poss*[numIters];
  for (int dupNum = 0; dupNum < numIters; ++dupNum) {
    Poss *newPoss = new Poss;
    newPosses[dupNum] = newPoss;
    map[dupNum] = new NodeMap;
    newPoss->Duplicate(poss, *(map[dupNum]), true, false);
    newPoss->PatchAfterDuplicate(*(map[dupNum]), true);
    bool trash;
    if (newPoss->RemoveLoops(&trash))
      LOG_FAIL("replacement for throw call");
  }
  
  Poss *rootPoss = newPosses[0];


  unsigned int numIn = poss->m_inTuns.size();
  //  cout << "orig numIn " << numIn << endl;
  for(unsigned int tunNum = 0; tunNum < numIn; ++tunNum) {
    //    cout << "***Working on tunNum " << tunNum << endl;
    LoopTunnel *possTunIn = (LoopTunnel*)(poss->m_inTuns[tunNum]);
    //    cout << possTunIn << endl;

    Tunnel *newTun = NULL;

    NodeVec newOutTuns;
    NodeVec outTuns;

    NodeConnVecIter childIter = possTunIn->m_children.begin();
    for(; childIter != possTunIn->m_children.end(); ++childIter) {
      Node *child = (*childIter)->m_n;
      if (child->IsTunnel(POSSTUNOUT)
	  && !FoundInNodeVec(outTuns, child)) {
	outTuns.push_back(child);
      }
    }
    if (outTuns.size() > 1)
      LOG_FAIL("replacement for throw call");

    if (!possTunIn->IsSplit()) {
      newTun = new Tunnel(POSSTUNIN);

      Node *prev = NULL;
      NodeVec inputToOutTuns;
      vector<ConnNum> inputNumToOutTuns;
      bool allInputsAreTunnel = true;
      
      if (!outTuns.empty()) {
	NodeVecIter outIter = outTuns.begin();
	for(; outIter != outTuns.end(); ++outIter) {
	  Tunnel *outTun = (Tunnel*)(*outIter);
	  if (!outTun->IsTunnel(POSSTUNOUT))
	    LOG_FAIL("replacement for throw call");
	  if (outTun->IsCombine())
	    LOG_FAIL("replacement for throw call");
	  
	  Node *input = outTun->Input(0);
	  inputToOutTuns.push_back(input);
	  if (input != possTunIn)
	    allInputsAreTunnel = false;

	  inputNumToOutTuns.push_back(outTun->InputConnNum(0));
	}
      }

      if (outTuns.empty() || allInputsAreTunnel) {
	for(int dupNum = 0; dupNum < numIters; ++dupNum) {
	  Poss *dup = newPosses[dupNum];
	  LoopTunnel *possTunIn = (LoopTunnel*)(dup->m_inTuns[tunNum]);
	  RedirectChildrenThatArentTuns(possTunIn, 0, newTun, 0, POSSTUNOUT);
	  possTunIn->RemoveAllChildren2Way();
	}
      }
      else {
	if (outTuns.size() > 1)
	  LOG_FAIL("replacement for throw call");
	prev = NULL;
	ConnNum inputNumToOutTun = inputNumToOutTuns[0];
	Node *inputToOutTun = inputToOutTuns[0];
	for(int dupNum = 0; dupNum < numIters; ++dupNum) {
	  Poss *dup = newPosses[dupNum];
	  LoopTunnel *possTunIn = (LoopTunnel*)(dup->m_inTuns[tunNum]);
	  if (dupNum) {
	    RedirectChildrenThatArentTuns(possTunIn, 0, prev, inputNumToOutTun, POSSTUNOUT);
	  }
	  else {
	    RedirectChildrenThatArentTuns(possTunIn, 0, newTun, 0, POSSTUNOUT);
	  }
	  //	  RedirectChildrenThatArentTuns(possTunIn, 0, newTun, 0);
	  NodeMapIter find = map[dupNum]->find(inputToOutTun);
	  if (find == map[dupNum]->end()) {
	    LOG_FAIL("replacement for throw call");
	  } else {
	    prev = find->second;
	  }
	}
      }
      
      if (!outTuns.empty()) {
	if (prev) {
	  if (outTuns.size() > 1)
	    LOG_FAIL("replacement for throw call");
	  Node *newOutTun = new Tunnel (POSSTUNOUT);
	  newOutTuns.push_back(newOutTun);
	  newOutTun->AddInput(prev, inputNumToOutTuns[0]);
	}
	else {
	  NodeVecIter outTunsIter = outTuns.begin();
	  for(; outTunsIter != outTuns.end(); ++outTunsIter) {
	    Node *newOutTun = new Tunnel (POSSTUNOUT);
	    newOutTuns.push_back(newOutTun);
	    newOutTun->AddInput(newTun, 0);
	  }
	}
      }
    }
    else { // tunnel is a split
      newTun = new Tunnel(POSSTUNIN);

      if (possTunIn->GetNodeClass() != SplitSingleIter::GetClass())
	LOG_FAIL("replacement for throw call");
      SplitSingleIter *split = (SplitSingleIter*)possTunIn;

      if (split->m_dir != PARTDOWN &&
	  split->m_dir != PARTRIGHT)
	LOG_FAIL("replacement for throw call");
      
      ViewMultipleIters *view = new ViewMultipleIters(split->m_dir,
						      split->GetMyLoop()->GetBSSize(),
						      numIters);

      CombineMultipleIters *com = new CombineMultipleIters(split->m_dir,
							   split->GetMyLoop()->GetBSSize(),
							   numIters);

      view->AddInput(newTun, 0);
      rootPoss->AddNode(view);
      rootPoss->AddNode(com);

      Node *inputToOutTun = NULL;
      ConnNum inputNumToOutTun = 99999;

      if (outTuns.empty())
	LOG_FAIL("replacement for throw call");
      if (outTuns.size() > 1)
	LOG_FAIL("replacement for throw call");
      Tunnel *outTun = (Tunnel*)(outTuns[0]);
      if (outTun->GetNodeClass() != CombineSingleIter::GetClass())
	LOG_FAIL("replacement for throw call");
      if (!outTun->IsTunnel(POSSTUNOUT))
	LOG_FAIL("replacement for throw call");

      inputToOutTun = outTun->Input(1);
      inputNumToOutTun = outTun->InputConnNum(1);

      bool passThrough = inputToOutTun->IsTunnel(POSSTUNIN);
      
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


	RedirectChildrenThatArentTuns(possTunIn, 1, view, dupNum, POSSTUNOUT);



	NodeConnVecIter iter = possTunIn->m_children.begin();
	for( ; iter != possTunIn->m_children.end(); ++iter) {
	  //part down and part right have 3 outputs for the split
	  // the 3rd (0-based) is the entire input
	  //this checks that all children only used the 1st (0-based)
	  // output 
	  if ((*iter)->m_num != 3) {
	    if (!(*iter)->m_n->IsTunnel(POSSTUNOUT)) {
	      cout << "num " << (*iter)->m_num << endl;
	      LOG_FAIL("replacement for throw call");
	    }
	  }
	}
	//This means it's a split and combine of read-only data
	// so it's just passed through
	//That's handled below
	if (!passThrough) {
	  NodeMapIter find = map[dupNum]->find(inputToOutTun);
	  if (find == map[dupNum]->end()) {
	    LOG_FAIL("replacement for throw call");
	  } else {
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

      Node *newOutTun = new Tunnel (POSSTUNOUT);
      newOutTuns.push_back(newOutTun);
      newOutTun->AddInput(com, 0);
    }

    rootPoss->AddNode(newTun);
    rootPoss->m_inTuns.push_back(newTun);

    if (!newOutTuns.empty()) {
      if (newOutTuns.size() != outTuns.size())
	LOG_FAIL("replacement for throw call");
      NodeVecIter newOutIter = newOutTuns.begin();
      NodeVecIter outIter = outTuns.begin();
      for(; outIter != outTuns.end(); ++outIter,++newOutIter) {
	Node *newOutTun = *newOutIter;
	Node *outTun = *outIter;
	rootPoss->AddNode(newOutTun);
	NodeVecIter outIter = rootPoss->m_outTuns.begin();
	Node *mappedOutTun = (*(map[0]))[outTun];
	if (!mappedOutTun)
	  LOG_FAIL("replacement for throw call");
	mappedOutTun->RemoveAllInputs2Way();
	for(; outIter != rootPoss->m_outTuns.end(); ++outIter) {
	  if (*outIter == mappedOutTun)
	    break;
	}
	if (outIter == rootPoss->m_outTuns.end())
	  LOG_FAIL("replacement for throw call");
	NodeVecIter inserted = rootPoss->m_outTuns.insert(outIter, newOutTun);
	inserted++;
	if (*inserted != mappedOutTun)
	  LOG_FAIL("replacement for throw call");
	rootPoss->m_outTuns.erase(inserted);
	rootPoss->DeleteNode(mappedOutTun);
	//      newOutTun->AddInput(newTun, 1);
	for(int dupNum = 1; dupNum < numIters; ++dupNum) {
	  NodeMapIter find = map[dupNum]->find(outTun);
	  if (find == map[dupNum]->end()) {
	    LOG_FAIL("replacement for throw call");
	  } else {
	    Node *possOutTun = find->second;
	    possOutTun->RemoveAllInputs2Way();
	    newPosses[dupNum]->DeleteNode(possOutTun);
	  }
	}
      }
    }
    else if (!outTuns.empty())
      LOG_FAIL("replacement for throw call");
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
	  if (!child->IsTunnel(POSSTUNOUT)) {
	    cout << poss->m_sets.size() << endl;
	    cout << child->GetNodeClass() << endl;
	    cout << child->GetType() << endl;
	    cout << "using output " << (*iter)->m_num << endl;
	    cout << "child " << child << endl;
	    cout << "tun " << tun << endl;
	    LOG_FAIL("replacement for throw call");
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
    Tunnel *outTun = (Tunnel*)(*outIter);
    if (outTun->IsLoopTunnel()) {
      cout << "found outTun " << outTun << endl;
      cout << outTun->GetNodeClass() << endl;
      cout << outTun->GetType() << endl;
      cout << "has " << outTun->m_inputs.size() << endl;
      for(int i = 0; i < outTun->m_inputs.size(); ++i) {
	cout << "input is " << outTun->m_inputs[i]->m_n << endl;
	cout << outTun->m_inputs[i]->m_n->GetType() << endl;
      }
      LOG_FAIL("replacement for throw call");
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
      if (node->IsTunnel()) {
	Tunnel *tun = (Tunnel*)node;
	if (tun->m_tunType == POSSTUNIN ||
	    tun->m_tunType == POSSTUNOUT)
	  LOG_FAIL("replacement for throw call");
      }
      node->m_poss = NULL;
      rootPoss->AddNode(node);
      NodeConnVecIter iter2 = node->m_children.begin();
      for(; iter2 != node->m_children.end(); ++iter2) {
	NodeConn *conn = *iter2;
	if (conn->m_n->IsLoopTunnel())
	  LOG_FAIL("replacement for throw call");
      }
    }
    dup->m_possNodes.clear();

    PSetVecIter iter2 = dup->m_sets.begin();
    for(; iter2 != dup->m_sets.end(); ++iter2) {
      BasePSet *set = *iter2;
      set->m_ownerPoss = NULL;
      rootPoss->AddPSet(set, true);
    }
    dup->m_sets.clear();
    
    delete dup;
  }

  if (rootPoss->m_outTuns.size() != poss->m_outTuns.size()) {
    cout << "bad size\n";
    LOG_FAIL("replacement for throw call");
  }

  delete [] newPosses;
  delete [] map;


  return rootPoss;
}

void FullyUnrollLoop::Apply(Node *node) const
{
  if (node->GetNodeClass() != SplitSingleIter::GetClass()) {
    cout << "Error: Unrolling node other than SplitSingleIter\n";
    LOG_FAIL("replacement for throw call");
  }
  
  SplitSingleIter *split = (SplitSingleIter*)node;
  
  LoopInterface *loop = split->GetMyLoop();
  BasePSet *base = dynamic_cast<BasePSet*>(loop);

  PossMMap &posses = base->GetPosses();

  PossMMapIter possIter = posses.begin();
  while (!possIter->second->ContainsNonLoopCode())
    ++possIter;

  Poss *rootPoss = UnrollPoss(possIter->second, loop, m_numIters);
  ++possIter;

  if (rootPoss->m_inTuns.size() != base->m_inTuns.size())
    LOG_FAIL("replacement for throw call");

  for (unsigned int i = 0; i < rootPoss->m_inTuns.size(); ++i) {
    Node *possInTun = rootPoss->m_inTuns[i];
    Node *loopInTun = base->m_inTuns[i];
    possInTun->AddInput(loopInTun->Input(0), loopInTun->InputConnNum(0));
  }

  RealPSet *newSet = new RealPSet(rootPoss);
  base->m_ownerPoss->AddPSet(newSet);

  if (newSet->m_outTuns.size() != base->m_outTuns.size()) {
    cout << newSet->m_outTuns.size() << endl;
    cout << base->m_outTuns.size() << endl;
    cout << rootPoss->m_outTuns.size() << endl;
    LOG_FAIL("replacement for throw call");
  }

  for(unsigned int i = 0; i < base->m_outTuns.size(); ++i) {
    Tunnel *loopTun = (Tunnel*)(base->m_outTuns[i]);
    if (!loopTun->m_children.empty()) {
      Tunnel *setTun = (Tunnel*)(newSet->m_outTuns[i]);
      loopTun->RedirectChildren(setTun);
    }
  }

  for(; possIter != posses.end(); ++possIter) {
    Poss *poss = possIter->second;
    if (poss->ContainsNonLoopCode()) {
      Poss *newPoss = UnrollPoss(poss, loop, m_numIters);
      newSet->AddPoss(newPoss);
    }
  }

  TunVecIter inIter = base->m_inTuns.begin();
  for( ; inIter != base->m_inTuns.end(); ++inIter) {
    Node *in = *inIter;
    in->RemoveAllInputs2Way();
    base->m_ownerPoss->DeleteNode(in);
  }
  base->m_inTuns.clear();



  for(unsigned int i = 0; i < base->m_outTuns.size(); ++i) {
    Tunnel *loopTun = (Tunnel*)(base->m_outTuns[i]);
    if (!loopTun->m_children.empty())
      LOG_FAIL("replacement for throw call");
    base->m_ownerPoss->DeleteNode(loopTun);
  }
  base->m_outTuns.clear();

  base->m_ownerPoss->RemoveFromSets(base);

  delete base;

  if (newSet->m_posses.empty()) {
    LOG_FAIL("replacement for throw call");
  }
}

CompactlyUnrollLoop::CompactlyUnrollLoop(int maxNumIters) 
  : m_numIters(maxNumIters) 
{
  if (m_numIters <= 0)
    LOG_FAIL("replacement for throw call");
}

bool CompactlyUnrollLoop::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != SplitSingleIter::GetClass()) {
    cout << "Error: Attempted to apply FullyUnrollLoop to non SplitSingleIter node\n";
    LOG_FAIL("replacement for throw call");
  }

  if (!node->IsTunnel(SETTUNIN))
    return false;

  const SplitSingleIter *split = (SplitSingleIter*)node;

  if (!split->m_isControlTun)
    return false;

  const BasePSet *loop = dynamic_cast<const BasePSet*>(split->GetMyLoop());
  if (!loop->IsReal()) {
    cout << "might be able to apply, but this is a shadow\n";
    LOG_FAIL("replacement for throw call");
  }

  if (loop->m_flags & SETLOOPISUNROLLED)
    return false;
  
  unsigned int numExecs = split->NumberOfLoopExecs();
  if (!numExecs) {
    cout << "Error: Attempted to unroll loop with 0 executions\n";
    LOG_FAIL("replacement for throw call");
  }

  for(unsigned int i = 0; i < numExecs; ++i) {
    unsigned int numIters = split->NumIters(i);
    if (numIters != m_numIters)
      return false;
  }

  const PossMMap &map = loop->GetPosses();
  PossMMapConstIter iter = map.begin();
  for(; iter != map.end(); ++iter) {
    const Poss *poss = iter->second;
    if (poss->ContainsNonLoopCode()) {
      //      static int count = 0;
      //      if (count)
      //	return false;
      //      else
      //	++count;
      return true;
    }
  }
  return false;
}

void CompactlyUnrollLoop::Apply(Node *node) const
{
  if (node->GetNodeClass() != SplitSingleIter::GetClass()) {
    cout << "Error: Unrolling node other than SplitSingleIter\n";
    LOG_FAIL("replacement for throw call");
  }
  
  SplitSingleIter *split = (SplitSingleIter*)node;
  
  LoopInterface *loop = split->GetMyLoop();
  BasePSet *base = dynamic_cast<BasePSet*>(loop);
  base->m_flags |= SETLOOPISUNROLLED;

  ((RealPSet*)base)->m_functionality += "unrolled\n";

  bool trash;
  ((RealPSet*)base)->RemoveLoops(&trash);
}

  
PartiallyUnrollLoop::PartiallyUnrollLoop(unsigned int unrollingFactor)
  : m_unrollingFactor(unrollingFactor)
{
}

bool PartiallyUnrollLoop::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != SplitSingleIter::GetClass()) {
    cout << "Error: Attempted to apply FullyUnrollLoop to non SplitSingleIter node\n";
    LOG_FAIL("replacement for throw call");
  }

  if (!node->IsTunnel(SETTUNIN))
    return false;

  const SplitSingleIter *split = (SplitSingleIter*)node;

  if (!split->m_isControlTun)
    return false;

  const LoopInterface *loopInt = split->GetMyLoop();
  const BasePSet *loop = dynamic_cast<const BasePSet*>(loopInt);
  if (!loop->IsReal()) {
    cout << "might be able to apply, but this is a shadow\n";
    LOG_FAIL("replacement for throw call");
  }

  if (loopInt->GetBSSize().m_multiple != 1)
    return false;

  if (loop->m_flags & SETLOOPISUNROLLED)
    return false;
  
  unsigned int numExecs = split->NumberOfLoopExecs();
  if (!numExecs) {
    cout << "Error: Attempted to unroll loop with 0 executions\n";
    LOG_FAIL("replacement for throw call");
  }

  //  unsigned int size = m_unrollingFactor * loopInt->GetBS();
  for(unsigned int i = 0; i < numExecs; ++i) {
    unsigned int numIters = split->NumIters(i);
    if (numIters % m_unrollingFactor)
      return false;
  }

  const PossMMap &map = loop->GetPosses();
  PossMMapConstIter iter = map.begin();
  for(; iter != map.end(); ++iter) {
    const Poss *poss = iter->second;
    if (poss->ContainsNonLoopCode()) {
      return true;

    }
  }
  return false;
}

void PartiallyUnrollLoop::Apply(Node *node) const
{
  if (node->GetNodeClass() != SplitSingleIter::GetClass()) {
    cout << "Error: Unrolling node other than SplitSingleIter\n";
    LOG_FAIL("replacement for throw call");
  }
  
  SplitSingleIter *split = (SplitSingleIter*)node;
  
  LoopInterface *loop = split->GetMyLoop();
  RealLoop *innerLoop = dynamic_cast<RealLoop*>(loop);

  Poss *outerPoss = innerLoop->m_ownerPoss;

  bool trash;
  innerLoop->RemoveLoops(&trash);
  
  outerPoss->RemoveFromSets(innerLoop);

  RealLoop *newOuter = (RealLoop*)(innerLoop->GetNewInst());
  newOuter->m_functionality = innerLoop->m_functionality + "partunrolled";

  TunVec newInnerPossInTuns;


  TunVec newSetTunsIn;
  TunVecIter tunIter = innerLoop->m_inTuns.begin();
  for(; tunIter != innerLoop->m_inTuns.end(); ++tunIter) {
    LoopTunnel *oldTun = (LoopTunnel*)(*tunIter);
    LoopTunnel *newTun = (LoopTunnel*)(oldTun->GetSetTunnel());
    newTun->m_tunType = POSSTUNIN;

    LoopTunnel *newSetTun = (LoopTunnel*)(oldTun->GetSetTunnel());


    newSetTun->m_pset = newOuter;
    newSetTun->m_tunType = SETTUNIN;
    newSetTun->AddInput(oldTun->Input(0), oldTun->InputConnNum(0));
    oldTun->RemoveAllInputs2Way();
    newSetTunsIn.push_back(newSetTun);
    newInnerPossInTuns.push_back(newTun);
    


    if (oldTun->IsSplit()) {
      if (oldTun->NumOutputs() != 4)
	LOG_FAIL("replacement for throw call");
      oldTun->AddInput(newTun, 1);
    }
    else
      oldTun->AddInput(newTun, 0);
    outerPoss->RemoveFromGraphNodes(oldTun);
  }


  TunVec newSetTunsOut;
  NodeVec newPossTunsOut;
  tunIter = innerLoop->m_outTuns.begin();
  for(; tunIter != innerLoop->m_outTuns.end(); ++tunIter) {
    LoopTunnel *oldTun = (LoopTunnel*)(*tunIter);
    LoopTunnel *newTun = (LoopTunnel*)(oldTun->GetSetTunnel());
    newTun->m_tunType = POSSTUNOUT;
    newPossTunsOut.push_back(newTun);
    LoopTunnel *newPossTunIn = (LoopTunnel*)(oldTun->GetMatchingInTun()->Input(0));


    LoopTunnel *newSetTun = (LoopTunnel*)(oldTun->GetSetTunnel());
    newSetTun->m_pset = newOuter;
    newSetTun->m_tunType = SETTUNOUT;
    //    newSetTun->AddInput(oldTun->Input(0), oldTun->InputConnNum(0));
    if (oldTun->m_children.size())
      oldTun->RedirectChildren(newSetTun);
    newSetTunsOut.push_back(newSetTun);
    
    if (newTun->IsCombine()) {
      newTun->AddInput(newPossTunIn, 0);
      newTun->AddInput(oldTun, 0);
      newTun->AddInput(newPossTunIn, 2);
      newTun->AddInput(newPossTunIn, 3);
    }
    else {
      newTun->AddInput(oldTun, 0);
      newTun->AddInput(newPossTunIn, 1);
    }
      
    outerPoss->RemoveFromGraphNodes(oldTun);
  }

  swap(newOuter->m_inTuns,newSetTunsIn);
  swap(newOuter->m_outTuns,newSetTunsOut);

  innerLoop->m_functionality += "innerUnrolled";

  Poss *newInnerPoss = new Poss(newPossTunsOut, true, true);

  if (newInnerPoss->m_inTuns.size() != newInnerPossInTuns.size())
    LOG_FAIL("replacement for throw call");

  newInnerPossInTuns.reserve(newInnerPoss->m_inTuns.size());
  for(auto tun : newInnerPossInTuns)
    newInnerPoss->m_inTuns.push_back(tun);
  newInnerPoss->m_inTuns.empty();
  //  swap(newInnerPoss->m_inTuns, newInnerPossInTuns);

  newOuter->AddPoss(newInnerPoss);
  outerPoss->AddPSet(newOuter, true, true);


  innerLoop->m_flags |= SETLOOPISUNROLLED;

#if TWOD
  newOuter->SetDimName(innerLoop->GetDimName());
#endif
  newOuter->m_bsSize = innerLoop->m_bsSize;
  newOuter->m_bsSize.m_multiple = m_unrollingFactor;

  newOuter->m_ownerPoss->ClearBeforeProp();
  newOuter->m_ownerPoss->ClearDataTypeCache();
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
      LOG_FAIL("replacement for throw call");
    Input(0)->Prop();
    if (m_numIters <= 1)
      LOG_FAIL("replacement for throw call");
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
      m_sizes->AddRepeatedSizes(m_bs.GetSize(), GetInputM(0)->NumSizes(), 1);
      m_info.m_numRowsVar = m_bs.VarName();
      break;
    }
  case (PARTRIGHT):
    {
      m_sizes->AddRepeatedSizes(m_bs.GetSize(), GetInputN(0)->NumSizes(), 1);
      m_info.m_numColsVar = m_bs.VarName();
      break;
    }
  default:
    LOG_FAIL("replacement for throw call");
  }
}

const Sizes* ViewMultipleIters::GetM(ConnNum num) const
{
  if (num > m_numIters) 
    LOG_FAIL("replacement for throw call");
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
    LOG_FAIL("replacement for throw call");
}

const Sizes* ViewMultipleIters::GetN(ConnNum num) const
{
  if (num > m_numIters)
    LOG_FAIL("replacement for throw call");

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
    LOG_FAIL("replacement for throw call");
}

Name ViewMultipleIters::GetName(ConnNum num) const
{
  if (num > m_numIters)
    LOG_FAIL("replacement for throw call");

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
      LOG_FAIL("replacement for throw call");
    }
  name.m_name += std::to_string((long long int) num);
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
      LOG_FAIL("replacement for throw call");
    *out << i << " = " << inName;


    if (i) {
      *out << " + " << std::to_string((long long int) i) << " * ";
      if (m_partDir == PARTDOWN) {
	if (!IsUnitStride(type.m_rowStride))
	  *out << type.m_rowStrideVar << " * ";
      }
      else if (m_partDir == PARTRIGHT) {
	if (!IsUnitStride(type.m_colStride))
	  *out << type.m_colStrideVar << " * ";
      }
      else
	LOG_FAIL("replacement for throw call");
      *out << m_bs.VarName();

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
    LOG_FAIL("replacement for throw call");
  for(int i = 0; i < m_numIters; ++i) {
    Var var(inName, i, m_info.m_type);
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


const DataTypeInfo& ViewMultipleIters::DataType(ConnNum num) const
{
  if (num < m_numIters) {
    if (!m_sizes)
      LOG_FAIL("replacement for throw call");
    return m_info;
  }
  else
    return InputDataType(0);
}


void CombineMultipleIters::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (m_inputs.size() != m_numIters+1)
      LOG_FAIL("replacement for throw call");

    if (m_numIters <= 1)
      LOG_FAIL("replacement for throw call");

    for (int i = 0; i < m_numIters+1; ++i)
      Input(i)->Prop();

    Node *in = Input(m_numIters);
    if (in->GetNodeClass() != ViewMultipleIters::GetClass())
      LOG_FAIL("replacement for throw call");
    ViewMultipleIters *inView = (ViewMultipleIters*)in;
    if (inView->m_bs != m_bs 
	|| inView->m_numIters != m_numIters
	|| inView->m_partDir != m_partDir)
      LOG_FAIL("replacement for throw call");


    switch (m_partDir) {
    case(PARTDOWN):
      {
	Size bs = m_bs.GetSize();
	for (int i = 0; i < m_numIters; ++i) {
	  if (*GetInputM(i) != bs) {
	    (*GetInputM(i)).Print();
	    LOG_FAIL("replacement for throw call");
	  }
	}	  
	break;
      }
    case(PARTRIGHT):
      {
	Size bs = m_bs.GetSize();
	for (int i = 0; i < m_numIters; ++i) {
	  if (*GetInputN(i) != bs) {
	    cout << "m_numIters " << m_numIters << endl;
	    for (int j = 0; j <= m_numIters; ++j) {
	      cout << "N for input j = " << j << endl;
	      (*GetInputN(j)).Print();
	    }
	    (*GetInputN(i)).Print();
	    cout << bs << endl;
	    LOG_FAIL("replacement for throw call");
	  }
	}
	break;
      }
    default:
      LOG_FAIL("replacement for throw call");
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


const Sizes* CombineMultipleIters::GetM(ConnNum num) const
{
  if (num > 0) 
    LOG_FAIL("replacement for throw call");
  return GetInputM(m_numIters);
}

const Sizes* CombineMultipleIters::GetN(ConnNum num) const
{
  if (num > 0) 
    LOG_FAIL("replacement for throw call");
  return GetInputN(m_numIters);
}

Name CombineMultipleIters::GetName(ConnNum num) const
{
  if (num > 0) 
    LOG_FAIL("replacement for throw call");
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

const DataTypeInfo& CombineMultipleIters::DataType(ConnNum num) const
{
  return InputDataType(m_numIters);
}

#endif


