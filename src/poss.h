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



#pragma once

#include "base.h"
#include "node.h"
#include "DLANode.h"
#include "basePSet.h"
#include "possTunnel.h"
#include <stdarg.h>
#include <unordered_set>
#include "realPSet.h"

//using namespace __gnu_cxx;

class Node;
class RealPSet;

class Poss
{
  static GraphNum M_count;
  size_t m_hash;
  bool m_hashValid;
 public:
  NodeVec m_possNodes;
  bool m_isSane;
  GraphNum m_parent;
  GraphNum m_num;
  NodeVec m_inTuns;
  NodeVec m_outTuns;
  RealPSet *m_pset;
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
  void PatchAfterDuplicate(NodeMap &map, bool deleteSetTunConnsIfMapNotFound = false);
  void DeleteChildAndCleanUp(Node *output, bool GoThroughTunnels=false, bool handleTunnelsAsNormalNodes=false, bool stopAtPossTunnels=false);
  virtual void DeleteNode(Node *node);
  void ForcePrint();
  bool CanPrint() const;
  virtual bool IsBoundary(Node *node) {return node->IsPossTunnel();}
  virtual void Duplicate(const Poss *orig, NodeMap &map, bool possMerging, bool useShadows);
  bool operator==(Poss &rhs);
  bool operator!=(Poss &rhs) {return !(*this == rhs);}
  virtual bool IsSane() const {return m_isSane;}
  void AddNode(Node *node);
  void TakeOverNode(Node *node);
  void AddNodes(int numNodes, ...);
  virtual void AddPSet(BasePSet *pset, bool expectToBeNew);
  virtual void AddUp(NodeVec &vec, Node *node, bool start, bool disconnectFromOwner);
  virtual void AddPSet(BasePSet *pset);
  bool ContainsNonLoopCode() const;
  bool RemoveLoops(bool *doneSomething);
  bool Simplify(const TransMap &simplifiers, bool recursive = false);
  void PrintTransVec();
  void RemoveConnectionToSet();
  void ExpandTunnels();
  Node* InTun(unsigned int num) const {return m_inTuns[num];}
  Node* OutTun(unsigned int num) const {return m_outTuns[num];}
  void MergePart1(unsigned int left, unsigned int right, 
		  BasePSet **leftSet, BasePSet **rightSet,
		  CullFunction cullFunc);
  void MergePart2(BasePSet *newSet, unsigned int left, NodeMap &map);
  void MergePart4(BasePSet *newSet, 
		  BasePSet *leftSet, 
		  BasePSet *rightSet, 
		  NodeMap &map,
		  NodeVec &newInputTunnelsToFix);
  void MergePart6(BasePSet *newSet, BasePSet *leftSet, 
		  BasePSet *rightSet, NodeMap &map);
  bool MergePosses(PossMMap &newPosses, const TransMap &simplifiers, CullFunction cullFunc);
  void MergePosses(unsigned int left, unsigned int right, const TransMap &simplifiers, CullFunction cullFunc);
  void FormSets(unsigned int phase);
  void FuseLoops(unsigned int left, unsigned int right, const TransMap &simplifiers, CullFunction cullFunc);
  virtual void Prop();
  virtual void Cull(Phase phase);
  virtual void ClearBeforeProp();
  virtual void ClearNodesPrinted();
  void ClearFullyExpanded();
  virtual void ClearPrintedFromGraph();
  string GetFunctionalityString() const;
  GraphNum TotalCount() const;
  bool TakeIter(const TransMap &transMap, const TransMap &simplifiers, 
		PossMMap &newPosses);
  bool HasFused(const BasePSet *left, const BasePSet *right) const;
  void SetFused(const BasePSet *left, const BasePSet *right);
  void RemoveFromGraphNodes(Node *node);
  void RemoveFromSets(BasePSet *set);
  void PrintNodeAddresses() const;
  BasePSet* FormSubPSet(NodeVec &outputTuns, bool isCritSect);
  void FillClique(NodeSet &set);
  BasePSet* FormSetForClique(NodeSet &set, bool isCritSect);

  static size_t Hash(const string &str);
  size_t GetHash();
  virtual void InvalidateHash() {m_hashValid=false;}

  void BuildDataTypeCache();
  void ClearDataTypeCache();

  void Flatten(ofstream &out) const;
  static void FlattenStatic(ofstream &out);
  void Unflatten(ifstream &in, SaveInfo &info);
  static void UnflattenStatic(ifstream &in);

  void PrintSetConnections();
};

void AddUsersOfLiveOutput(Node *node, ConnNum connNum, NodeSet &set);
