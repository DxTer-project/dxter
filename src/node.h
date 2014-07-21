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
#define BUILDFLAG (1L<<1)
#define PRINTEDFLAG (1L<<2)
#define HASREFINEDFLAG (1L<<3)

class DataTypeInfo;
class Loop;

class Node
{
 public:
  Flags m_flags;
  TransSet m_applications; //Don't duplicate!
  TransSet m_inverseOps; //Duplicate!

  NodeConnVec m_inputs;
  NodeConnVec m_children;
  Poss *m_poss;

#ifdef TRACKORIG
  const Node *m_orig;
#endif

  //Implement at least these in sub-classes
  /*****************/
  //Get a new instance of the node's class 
  // (with default values or whatever since Duplicate will be called)
  virtual Node* GetNewInst() = 0;
  //Propagate data type information, calculate node costs, and check for
  // error conditions
  virtual void Prop() = 0;
  //The maximum phase in which this node should be found
  // Going from that phase to the next, a Poss with this node
  // should be culled (but HasRefined will be checked
  // and an exception will be thrown if the node hasn't refined)
  //Should super message
  virtual Phase MaxPhase() const {return NUMPHASES;}
  //Print code
  virtual void PrintCode(IndStream &out) = 0;
  //Duplicate all of the original node's type information
  // into me
  //Should super message
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  //Get the node type used for comparing two nodes for equality
  virtual NodeType GetType() const;
  //Number of outputs expected
  virtual unsigned int NumOutputs() const {return 1;}
  //Cost of the node (guaranteed to be called after Prop is called
  virtual Cost GetCost() = 0;
  //Get the node's class so you know its type
  virtual ClassType GetNodeClass() const = 0;
  //Static class type
  static ClassType GetClass() {return "node";}
  //Get the name of the output variable
  virtual Name GetName(unsigned int num) const = 0;
  //Returns true if the node overwrites the variable passed in by Node input
  virtual bool Overwrites(const Node *input, unsigned int num) const = 0;
  //Returns true if the variable passed in by Node input is also passed
  // out of the node
  virtual bool KeepsInputVarLive(Node *input, unsigned int numIn, 
				 unsigned int &numOut) const = 0;
  //Add any variable declarations for this node (e.g., new
  // variables that are used as temporaries)
  //Should super message
  virtual void AddVariables(VarSet &set) const;
  //Clear any cached information before propagating (e.g., cost
  // or data type info)
  virtual void ClearBeforeProp() {}
  //Clear the size cache specifically
  virtual void ClearDataTypeCache() {}
  //Build the size cache
  virtual void BuildDataTypeCache() {}
  /*****************/




  Node();
  virtual ~Node();
  void Cull(Phase phase);
  void CheckConnections();
  void AddChild(Node *node, unsigned int num);
  void RemoveChild(Node *node, unsigned int num);
  void RemoveInput(Node *node, unsigned int num);
  void RemoveAllInputs2Way();
  void RemoveAllChildren2Way();
  virtual bool operator==(const Node &rhs) const;
  virtual bool operator!=(const Node &rhs) const {return !(*this == rhs);}
  void PatchAfterDuplicate(NodeMap &map, bool deleteSetTunConnsIfMapNotFound = false);
  virtual void Print(IndStream &out, GraphNum graphNum);
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
  inline void ClearPrinted() {m_flags &= ~PRINTEDFLAG;}
  inline void SetPrinted() {m_flags |= PRINTEDFLAG;} 
  inline bool HasPrinted() const {return m_flags & PRINTEDFLAG;}
  inline void ClearHasRefined() {m_flags &= ~HASREFINEDFLAG;}
  inline void SetHasRefined() {m_flags |= HASREFINEDFLAG;} 
  inline bool HasRefined() const {return m_flags & HASREFINEDFLAG;}
  Node* Input(unsigned int num) const;
  NodeConn* InputConn(unsigned int num) const;
  //output number of the input taken as input
  unsigned int InputConnNum(unsigned int num) const;
  Node* Child(unsigned int num) const;
  //child of my output number num
  unsigned int ChildConnNum(unsigned int num) const;
  void AddToPoss(Poss *poss);
  unsigned int NumChildrenOfOutput(unsigned int num) const;
  bool InChildren(Node *node, unsigned int num) const;
  bool InInputs(Node *node, unsigned int num) const;
  virtual bool IsDLA() const {return false;}
  virtual bool IsPossTunnel() const {return false;}
  virtual bool IsPossTunnel(PossTunType type) const {return false;}
  virtual bool IsLoopTunnel() const {return false;}
  virtual bool IsParallel() const {return false;}
  const Loop* FindClosestLoop() const;

#if DOBLIS
  virtual Comm ParallelComm() const {throw;}
  virtual Comm WithinParallelism() const;
  virtual Comm HasBarrier() const {return CORECOMM;}
  virtual bool RemoveParallelization() {throw;}
  bool InCriticalSection() const;
#endif

  virtual bool IsDataDependencyOfInput() const {return true;}

  void PrintChildren();
  void PrintInputs();

  void BuildDataTypeCacheRecursive();

  void Flatten(ofstream &out) const;
  virtual void FlattenCore(ofstream &out) const = 0;
  void Unflatten(ifstream &in, SaveInfo &info);
  virtual void UnflattenCore(ifstream &in, SaveInfo &info) = 0;

  virtual const DataTypeInfo& DataType(unsigned int num) const = 0;
  virtual const DataTypeInfo& InputDataType(unsigned int num) const;

  string GetFunctionalityString() const;
};

void FullyFlatten(const NodeVec &vec, ofstream &out);
void FullyUnflatten(NodeVec &vec, ifstream &in, SaveInfo &info);
