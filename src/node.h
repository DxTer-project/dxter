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
#include "var.h"
#include <stdarg.h>
#include "comm.h"
#include "sizes.h"

#define NODIST
#define NOALGS

#define PRINTEMPTY 1


void PrintSetOrNodeChildren(Node *node);
void PrintSetOrNodeInputs(Node *node);


typedef unsigned int ConnNum;

//Keeps track of child/parent node and the
// input/output number
class NodeConn {
 public:
  Node *m_n;
  ConnNum m_num;
 NodeConn() : m_n(NULL) {}
  NodeConn(Node *n, ConnNum num);
  NodeConn(const NodeConn *conn);
  bool operator==(const NodeConn &rhs) const;
  void SetNode(Node *node);
  void Flatten(ofstream &out) const;
  void Unflatten(ifstream &in);
};

#define NODEBUILDFLAG (1L<<1)
#define NODEHASREFINEDFLAG (1L<<2)

class DataTypeInfo;
class RealLoop;
class ShadowLoop;
class BasePSet;
class GraphIter;

class Node
{
 public:
  Flags m_flags;
  TransSet m_applications; //Don't duplicate!
  TransSet m_inverseOps; //Duplicate!

  NodeConnVec m_inputs;
  NodeConnVec m_children;
  Poss *m_poss;

  //Implement at least these in subclasses
  /*****************/
  //Get a new instance of the node's class 
  // (with uninitialized values)
  virtual Node* GetNewInst() = 0;
  //Propagate data type information, calculate node costs, and check for
  // error conditions
  //Should super message
  virtual void Prop() = 0;
  //The maximum phase in which this node should be found
  // Going from that phase to the next, a Poss with this node
  // should be culled (but HasRefined will be checked
  // and an exception will be thrown if the node hasn't been refined)
  virtual Phase MaxPhase() const {return NUMPHASES;}
  //Print code
  virtual void PrintCode(IndStream &out) = 0;
  //Duplicate all of the original node's type information
  // into me
  //Should super message
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  //Get the node type used for comparing two nodes for equality
  // This should include enough information to compare two nodes
  // of the same class that have different parameters like
  // coefficients or transposition
  virtual NodeType GetType() const;
  //Number of outputs expected
  virtual unsigned int NumOutputs() const {return 1;}
  //Cost of the node (guaranteed to be called after Prop and BuildDataTypeCache is called
  virtual Cost GetCost() = 0;
  //Get the node's class so you know its type and can cast it
  // correctly from Node* to the most specific class
  virtual ClassType GetNodeClass() const = 0;
  //Static class type
  //Enables a comparison like 
  //  if (node->GetNodeClass() == FooClass::GetClass()) {
  //     FooClass *foo = (FooClass*)node;
  //     ...  do foo stuff  ...
  //  }
  static ClassType GetClass() {return "node";}
  //Get the name of the output variable
  virtual Name GetName(ConnNum num) const = 0;
  //Returns true if the node overwrites the variable passed in by Node input
  //This enables a dataflow graph, which shouldn't overwrite inputs,
  // to be efficiently output in code by reusing memory of inputs
  virtual bool Overwrites(const Node *input, ConnNum num) const = 0;
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
  //Size of the output
  virtual int Outputs() { return 0; }

  virtual const DataTypeInfo& DataType(ConnNum num) const = 0;

  /*****************/

  Node();
  virtual ~Node();
  void Cull(Phase phase);
  void CheckConnections();
  void AddChild(Node *node, ConnNum num);
  void RemoveChild(Node *node, ConnNum num);
  void RemoveInput(Node *node, ConnNum num);
  void RemoveAllInputs2Way();
  void RemoveAllChildren2Way();
  virtual bool operator==(const Node &rhs) const;
  virtual bool operator!=(const Node &rhs) const {return !(*this == rhs);}
  void PatchAfterDuplicate(NodeMap &map, bool deleteSetTunConnsIfMapNotFound = false);
  void Print(IndStream &out);
  void AddInput(Node *node);
  void AddInputs(int numArgs, ...);
  void AddInputs0(int numArgs, ...);

  //num is the output number of the input that this is taking as input
  void AddInput(Node *node, ConnNum num);
  void ChangeInput1Way(Node *oldInput, ConnNum oldNum, Node *newInput, ConnNum newNum);
  void ChangeInput2Way(Node *oldInput, ConnNum oldNum, Node *newInput, ConnNum newNum);
  void RedirectChildren(Node *newInput);
  void RedirectChildren(Node *newInput, ConnNum newNum);
  void RedirectChildren(ConnNum oldNum, Node *newInput, ConnNum newNum);
  void RedirectAllChildren(Node *newInput);
  void RedirectChild(unsigned int childNum, Node *newInput, ConnNum newNum);
  bool HasApplied(const Transformation*) const;
  bool Applied(const Transformation*);
  void SetPoss(Poss *poss) {m_poss=poss;}
  string GetNameStr(ConnNum num) const {return GetName(num).str();}
  Name GetInputName(ConnNum num) const;
  string GetInputNameStr(ConnNum num) const {return GetInputName(num).str();}

  inline void ClearHasRefined() {m_flags &= ~NODEHASREFINEDFLAG;}
  inline void SetHasRefined() {m_flags |= NODEHASREFINEDFLAG;}
  inline bool HasRefined() const {return m_flags & NODEHASREFINEDFLAG;}
  Node* Input(ConnNum num) const;
  NodeConn* InputConn(ConnNum num) const;
  //output number of the input taken as input
  ConnNum InputConnNum(ConnNum num) const;
  Node* Child(unsigned int num) const;
  //child of my output number num
  ConnNum ChildConnNum(ConnNum num) const;
  void AddToPoss(Poss *poss);
  unsigned int NumChildrenOfOutput(ConnNum num) const;
  bool InChildren(Node *node, ConnNum num) const;
  bool InInputs(Node *node, ConnNum num) const;
  virtual bool IsDLA() const {return false;}
  virtual bool IsTunnel() const {return false;}
  virtual bool IsTunnel(TunType type) const {return false;}
  virtual bool IsLoopTunnel() const {return false;}
  virtual bool IsParallel() const {return false;}
  const BasePSet* FindClosestLoop() const;

#if DORQO
  virtual bool IsSortable() const {return false;}
  virtual bool IsJoin() const {return false;}
  virtual bool IsProjection() const {return false;}
#endif

  virtual Node* ChildOfOutputWithClass(ConnNum num, string nodeClass) const;
  virtual bool OutputHasChildOfClass(ConnNum num, string nodeClass) const;

#if DOBLIS
  virtual Comm ParallelComm() const {LOG_FAIL("replacement for throw call");}
  virtual Comm WithinParallelism() const;
  virtual Comm HasBarrier() const {return CORECOMM;}
  virtual bool RemoveParallelization() {LOG_FAIL("replacement for throw call");}
  bool InCriticalSection() const;
#endif

  virtual bool IsDataDependencyOfInput() const {return true;}

#if DOTENSORS
  virtual bool CreatesNewVars() const {return true;}
#endif

  void PrintChildren() const;
  void PrintInputs() const;

  void PrintChildren(IndStream& out) const;
  void PrintInputs(IndStream& out) const;

  void BuildDataTypeCacheRecursive();

  void Flatten(ofstream &out) const;
  virtual void FlattenCore(ofstream &out) const {throw;}
  void Unflatten(ifstream &in, SaveInfo &info);
  virtual void UnflattenCore(ifstream &in, SaveInfo &info) {throw;}
  
  virtual const DataTypeInfo& InputDataType(ConnNum num) const;

#if DOLLDLA
  virtual const Type GetDataType() const;
  virtual const int GetVecRegWidth() const;
#endif

  string GetFunctionalityString() const;
};

void FullyFlatten(const NodeVec &vec, ofstream &out);
void FullyUnflatten(NodeVec &vec, ifstream &in, SaveInfo &info);
