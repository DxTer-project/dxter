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
#include "var.h"
#include <stdarg.h>
#include "comm.h"

#define NODIST
#define NOALGS

void PrintSetOrNodeChildren(Node *node);
void PrintSetOrNodeInputs(Node *node);

//#define TRACKORIG

//Keeps track of child/parent node and the
// input/output number
class NodeConn {
 public:
  Node *m_n;
  unsigned int m_num;
 NodeConn() : m_n(NULL) {}
  NodeConn(Node *n, unsigned int num);
  NodeConn(const NodeConn *conn);
  bool operator==(const NodeConn &rhs) const;
  void SetNode(Node *node);
  void Flatten(ofstream &out) const;
  void Unflatten(ifstream &in);
};

typedef unsigned int Flags;
#define BUILDFLAG 1L

class Node
{
 private:
  unsigned int m_num;
  bool m_hasPrinted;

 public:
  Flags m_flags;
  TransSet m_applications; //Don't duplicate!
  TransSet m_inverseOps; //Duplicate!

  NodeConnVec m_inputs;
  NodeConnVec m_children;
  Poss *m_poss;
  bool m_hasRefined;

#ifdef TRACKORIG
  const Node *m_orig;
#endif

  //Implement these in sub-classes
  virtual Node* GetNewInst() = 0;
  virtual void Prop() = 0;
  virtual Phase MaxPhase() const {return NUMPHASES;}
  virtual bool HasProped() const = 0;
  virtual void PrintCode(IndStream &out) = 0;
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const;
  virtual void SanityCheck();
  virtual unsigned int NumOutputs() const {return 1;}
  virtual string GetCostStr() = 0;
  virtual Cost GetCost() = 0;
  virtual ClassType GetNodeClass() const = 0;
  static ClassType GetClass() {return "node";}
  virtual Name GetName(unsigned int num) const = 0;
  virtual bool Overwrites(const Node *input, unsigned int num) const = 0;
  virtual bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const = 0;
#if DOTENSORS
  virtual void AddVariables(VarSet &set) const;
#endif

  virtual void ClearSizeCache() {}
  virtual void BuildSizeCache() {}
  void BuildSizeCacheRecursive();

  Node();
  virtual ~Node();
  void Cull(Phase phase);
  void AddChild(Node *node, unsigned int num);
  void RemoveChild(Node *node, unsigned int num);
  void RemoveInput(Node *node, unsigned int num);
  virtual bool operator==(const Node &rhs) const;
  virtual bool operator!=(const Node &rhs) const {return !(*this == rhs);}
  void PatchAfterDuplicate(NodeMap &map);
  virtual void Print(IndStream &out, unsigned int graphNum);
  bool CanPrintCode() const;
  void AddInput(Node *node);
  void AddInputs(int numArgs, ...);
  //num is the output number of the input that this is taking as input
  void AddInput(Node *node, unsigned int num);
  void ChangeInput1Way(Node *oldInput, unsigned int oldNum, Node *newInput, unsigned int newNum);
  void ChangeInput2Way(Node *oldInput, unsigned int oldNum, Node *newInput, unsigned int newNum);
  void RedirectChildren(Node *newInput);
  void RedirectChildren(Node *newInput, unsigned int newNum);
  void RedirectChildren(unsigned int oldNum, Node *newInput, unsigned int newNum);
  void RedirectAllChildren(Node *newInput);
  void RedirectChild(unsigned int childNum, Node *newInput, unsigned int newNum);
  bool HasApplied(const Transformation*) const;
  bool Applied(const Transformation*);
  void SetPoss(Poss *poss) {m_poss=poss;}
  string GetNameStr(unsigned int num) const {return GetName(num).str();}
  Name GetInputName(unsigned int num) const;
  string GetInputNameStr(unsigned int num) const {return GetInputName(num).str();}
  inline void ClearPrinted() {m_hasPrinted=false;}
  inline void SetPrinted() {m_hasPrinted=true;}
  inline bool HasPrinted() const {return m_hasPrinted;}
  Node* Input(unsigned int num) const;
  NodeConn* InputConn(unsigned int num) const;
  //output number of the input taken as input
  unsigned int InputConnNum(unsigned int num) const;
  Node* Child(unsigned int num) const;
  //child of my output number num
  unsigned int ChildConnNum(unsigned int num) const;
  virtual void ClearBeforeProp() {}
  void AddToPoss(Poss *poss);
  unsigned int NumChildrenOfOutput(unsigned int num) const;
  bool InChildren(Node *node, unsigned int num) const;
  bool InInputs(Node *node, unsigned int num) const;
  virtual bool IsDLA() const {return false;}
  virtual bool IsPossTunnel() const {return false;}
  virtual bool IsPossTunnel(PossTunType type) const;
  virtual bool IsLoopTunnel() const {return false;}
  virtual bool IsParallel() const {return false;}
  virtual Comm ParallelComm() const {throw;}
  virtual Comm WithinParallelism() const;
  virtual Comm HasBarrier() const {return CORECOMM;}
  virtual bool RemoveParallelization() {throw;}
  bool InCriticalSection() const;
  virtual bool IsDataDependencyOfInput() const {return true;}

  void PrintChildren();
  void PrintInputs();

  void Flatten(ofstream &out) const;
  virtual void FlattenCore(ofstream &out) const = 0;
  void Unflatten(ifstream &in, SaveInfo &info);
  virtual void UnflattenCore(ifstream &in, SaveInfo &info) = 0;
};

void FullyFlatten(const NodeVec &vec, ofstream &out);
void FullyUnflatten(NodeVec &vec, ifstream &in, SaveInfo &info);
