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



#pragma once

#include "base.h"
#include "node.h"
#include "DLANode.h"
//#include "basePSet.h"
#include "tunnel.h"
#include <stdarg.h>
#include <unordered_set>
//#include "realPSet.h"
//#include "basePSet.h"

//using namespace __gnu_cxx;

class Node;
class RealPSet;
class BasePSet;

#define POSSISSANEFLAG (1L<<1)

typedef vector<NodeConn*> NodeConnVec;
typedef NodeConnVec::iterator NodeConnVecIter;
typedef NodeConnVec::const_iterator NodeConnVecConstIter;

int FindInNodeConnVec(const Node *node,
		      const ConnNum num,
		      const NodeConnVec &vec);


//When merging and keeping track of where tunnels get merged,
// we have to do some bookeeping to note what the inputs
// of tunnels are and the children of tunnes
//If two different tunnels are input to the same child,
// we need to look at the input number for that child
// as a differentiator
//For two inputs to tunnels, we don't care
typedef struct { 
  NodeConn m_conn;
  ConnNum m_num;
} NodeConnAndNum;

struct NodeConnAndNumComp {
  bool operator() (const NodeConnAndNum &lhs, const NodeConnAndNum &rhs) const
  {return lhs.m_conn.m_n < rhs.m_conn.m_n || (lhs.m_conn.m_n == rhs.m_conn.m_n && lhs.m_conn.m_num < rhs.m_conn.m_num)
      || (lhs.m_conn.m_n == rhs.m_conn.m_n && lhs.m_conn.m_num == rhs.m_conn.m_num && lhs.m_num < rhs.m_num); }
};

    

typedef map<NodeConnAndNum,vector<int>,NodeConnAndNumComp> NodeConnAndNumIntMap;
typedef NodeConnAndNumIntMap::iterator NodeConnAndNumIntMapIter;



class Poss
{
  static GraphNum M_count;
  size_t m_hash;
  bool m_hashValid;
 public:
  NodeVec m_possNodes;
  PSetVec m_sets;
  Flags m_flags;
  GraphNum m_parent;
  GraphNum m_num;
  TunVec m_inTuns;
  TunVec m_outTuns;
  RealPSet *m_pset;
  TransVec m_transVec;
  bool m_fullyExpanded;
  static StrSet M_fusedSets;
  Cost m_cost;
  Poss();
  virtual ~Poss();
  Poss(Tunnel *tun);
  Poss(Node *node, bool goUp=false);
  Poss(int numArgs, ...);
  Poss(const NodeVec &nodes, bool outTuns, bool disconnectFromOwner);
  void InitHelper(const NodeVec &nodes, bool outTuns, bool disconnectFromOwner);
  void MarkInsane(bool wrongPhase = false);
  void PatchAfterDuplicate(NodeMap &map, bool deleteSetTunConnsIfMapNotFound = false);
  void DeleteChildAndCleanUp(Node *output, bool GoThroughTunnels=false, bool handleTunnelsAsNormalNodes=false, bool stopAtTunnels=false);
  virtual void DeleteNode(Node *node);
  void ForcePrint();
  virtual bool IsBoundary(Node *node) {return node->IsTunnel();}
  virtual void Duplicate(const Poss *orig, NodeMap &map, bool possMerging, bool useShadows);
  bool operator==(Poss &rhs);
  bool operator!=(Poss &rhs) {return !(*this == rhs);}
  inline bool IsSane() const {return m_flags & POSSISSANEFLAG; }
  void AddNode(Node *node);
  void RecursiveAddNodes(Node* node);
  void TakeOverNode(Node *node);
  void AddNodes(int numNodes, ...);
  virtual void AddPSet(BasePSet *pset, bool expectToBeNew, bool addTuns = false);
  virtual void AddUp(NodeVec &vec, Node *node, bool start, bool disconnectFromOwner);
  virtual void AddPSet(BasePSet *pset);
  bool ContainsNonLoopCode() const;
  bool ContainsLoops() const;
  bool RemoveLoops(bool *doneSomething);
  bool Simplify(const TransMap &simplifiers, bool recursive = false);
  void PrintTransVec();
  void PrintTransVecUp();
  void RemoveConnectionToSet();
  //  void ExpandTunnels();
  Node* InTun(unsigned int num) const {return m_inTuns[num];}
  Node* OutTun(unsigned int num) const {return m_outTuns[num];}
  bool MergePart1(unsigned int left, unsigned int right, 
		  BasePSet **leftSet, BasePSet **rightSet);
  void MergePart2(RealPSet *newSet, 
		  BasePSet *leftSet, BasePSet *rightSet,
		  unsigned int left, NodeMap &mapLeft, NodeMap &mapRight,
		  NodeConnAndNumIntMap &inMap, NodeConnAndNumIntMap &outMap);
  void MergePart4(RealPSet *newSet, 
		  BasePSet *leftSet, 
		  BasePSet *rightSet, 
		  NodeMap &mapLeft, NodeMap &mapRight,
		  TunVec &newInputTunnelsToFix);
  void MergePart6(RealPSet *newSet, BasePSet *leftSet, 
		  BasePSet *rightSet, NodeMap &mapLeft, NodeMap &mapRight);
  void MergePart7(RealPSet *newSet, 
		  unsigned int numLeftInTuns, unsigned int numRightInTuns,
		  unsigned int numLeftOutTuns, unsigned int numRightOutTuns,
		  NodeConnAndNumIntMap &inMap, NodeConnAndNumIntMap &outMap);
  bool MergePosses(PossMMap &newPosses, const TransMap &simplifiers, CullFunction cullFunc);
  void MergePosses(unsigned int left, unsigned int right, const TransMap &simplifiers, CullFunction cullFunc);
  void FormSets(unsigned int phase);
  void FuseLoops(unsigned int left, unsigned int right, const TransMap &simplifiers, CullFunction cullFunc);
  virtual Cost Prop();
  virtual void Cull(Phase phase);
  void CullWorstPerformers(double percentToCull, int ignoreThreshold);
  void CullAllBut(int num);
  virtual void ClearBeforeProp();
  void ClearFullyExpanded();
  string GetFunctionalityString() const;
  GraphNum TotalCount() const;
  bool TakeIter(const TransMap &transMap, const TransMap &simplifiers, 
		PossMMap &newPosses);
  bool HasFused(const BasePSet *left, const BasePSet *right) const;
  void SetFused(const BasePSet *left, const BasePSet *right);
  void RemoveFromGraphNodes(Node *node);
  void RemoveFromSets(BasePSet *set);
  void PrintNodeAddresses() const;
  RealPSet* FormSubPSet(NodeVec &outputTuns, bool isCritSect);
  void FillClique(NodeSet &set);
  static size_t Hash(const string &str);
  size_t GetHash();
  virtual void InvalidateHash() {m_hashValid=false;}

  void BuildDataTypeCache();
  void ClearDataTypeCache();

  virtual void InlineAllSets();

  void Flatten(ofstream &out) const;
  static void FlattenStatic(ofstream &out);
  void Unflatten(ifstream &in, SaveInfo &info);
  static void UnflattenStatic(ifstream &in);

  void PrintSetConnections();
  void ReplaceShadowSetWithReal(unsigned int i);
  void RemoveAndDeleteNodes(NodeVec &vec);
  void RemoveAndDeleteNodes(TunVec &vec);
  void SetDeletingRecursively();
  void ClearDeletingRecursively();

#if DOTENSORS
  bool CheckDistTypeSet(StrSet &set) const;
  void GetDistTypeSet(StrSet &set) const;
  bool HasRedist() const;
#endif // DOTENSORS
};

void AddUsersOfLiveOutput(Node *node, ConnNum connNum, NodeSet &set);
