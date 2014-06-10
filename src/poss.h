/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2013, The University of Texas and Bryan Marker

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



#pragma once

#include "base.h"
#include "node.h"
#include "DLANode.h"
#include "pset.h"
#include "possTunnel.h"
#include <stdarg.h>
#include <unordered_set>

//using namespace __gnu_cxx;

class Node;

class Loop;

class Poss
{
  static unsigned int M_count;
  size_t m_hash;
  bool m_hashValid;
 public:
  bool m_hasPrinted;
  NodeVec m_possNodes;
  bool m_isSane;
  unsigned int m_parent;
  unsigned int m_num;
  NodeVec m_inTuns;
  NodeVec m_outTuns;
  PSet *m_pset;
  TransVec m_transVec;
  bool m_fullyExpanded;
  PSetVec m_sets;
  static StrSet M_fusedSets;
  Poss();
  virtual ~Poss();
  Poss(PossTunnel *tun);
  Poss(Node *node, bool goUp=false);
  Poss(int numArgs, ...);
  Poss(const NodeVec &nodes, bool outTuns, bool disconnectFromOwner);
  void InitHelper(const NodeVec &nodes, bool outTuns, bool disconnectFromOwner);
  void MarkInsane(bool wrongPhase = false);
  void PatchAfterDuplicate(NodeMap &map);
  void DeleteChildAndCleanUp(Node *output, bool GoThroughTunnels=false, bool handleTunnelsAsNormalNodes=false);
  virtual void DeleteNode(Node *node);
  virtual Cost EvalCurr(TransConstVec &transList);
  virtual Cost EvalAndSetBest();
  virtual void Print(IndStream &out, unsigned int &graphNum);
  virtual void EvalRoot(IndStream &out, unsigned int &graphNum, unsigned int whichGraph, unsigned int &optGraph, Cost &optCost);
  virtual void PrintRoot(IndStream &out, unsigned int &graphNum, unsigned int whichGraph);
  virtual void PrintCurrRoot(IndStream &out, const VarSet &set);
  void ForcePrint();
  bool CanPrint() const;
  virtual bool IsBoundary(Node *node) {return node->IsPossTunnel();}
  virtual void Duplicate(const Poss *orig, NodeMap &map, bool possMerging);
  bool operator==(Poss &rhs);
  bool operator!=(Poss &rhs) {return !(*this == rhs);}
  virtual bool IsSane() const {return m_isSane;}
  void AddNode(Node *node);
  void AddNodes(int numNodes, ...);
  virtual void AddPSet(PSet *pset, bool expectToBeNew);
  virtual void AddUp(NodeVec &vec, Node *node, bool start, bool disconnectFromOwner);
  virtual void AddLoop(Loop *loop);
  bool Simplify(const TransMap &simplifiers);
  void PrintTransVec();
  void RemoveConnectionToSet();
  void ExpandTunnels();
  Node* InTun(unsigned int num) const {return m_inTuns[num];}
  Node* OutTun(unsigned int num) const {return m_outTuns[num];}
  bool MergePosses(PossMMap &newPosses, const TransMap &simplifiers, CullFunction cullFunc);
  void MergePosses(unsigned int left, unsigned int right, const TransMap &simplifiers, CullFunction cullFunc);
  void FormSets(unsigned int phase);
  void FuseLoops(unsigned int left, unsigned int right, const TransMap &simplifiers, CullFunction cullFunc);
  virtual void Prop();
  virtual void Cull(Phase phase);
  virtual void ClearBeforeProp();
  virtual void ClearNodesPrinted();
  void ClearFullyExpanded();
  virtual void ClearPrinted();
  void GetTransVec(TransVec &transVec) const;
  void GetCurrTransVec(TransVec &transVec) const;
  unsigned int TotalCount() const;
  bool TakeIter(const TransMap &transMap, const TransMap &simplifiers, 
		PossMMap &newPosses);
  bool GlobalSimplification(const TransMap &globalSimplifiers, const TransMap &simplifiers);
  bool HasFused(const Loop *left, const Loop *right) const;
  void SetFused(const Loop *left, const Loop *right);
  void RemoveFromGraphNodes(Node *node);
  void RemoveFromSets(PSet *set);
  void PrintNodeAddresses() const;
  PSet* FormSubPSet(NodeVec &outputTuns, bool isCritSect);
  void FillClique(NodeSet &set);
  PSet* FormSetForClique(NodeSet &set, bool isCritSect);

#if DOTENSORS
  void AddCurrPossVars(VarSet &set) const;
#endif

  void ClearCurrPoss();
  bool IncrementCurrPoss();

  size_t GetHash();
  virtual void InvalidateHash() {m_hashValid=false;}

  void BuildSizeCache();
  void ClearSizeCache();

  void Flatten(ofstream &out) const;
  static void FlattenStatic(ofstream &out);
  void Unflatten(ifstream &in, SaveInfo &info);
  static void UnflattenStatic(ifstream &in);

  void PrintSetConnections();
};

void AddUsersOfLiveOutput(Node *node, unsigned int connNum, NodeSet &set);
