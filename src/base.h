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

#include "layers.h"
#include "sizes.h"

//Maximum number of refinement to use
// in a MultiTrans (set to something large
// if you don't want to use this heuristic)
#define MAXNUMBEROFREFINEMENTS 2

//Output cost code for Matlab vs. R
#define MATLAB

extern char END;
extern char START;

#define WRITE(var) out.write((char*)&var,sizeof(var))
#define READ(var) in.read((char*)&var,sizeof(var))

typedef unsigned int Phase;


inline bool IsValidCost(Cost cost)
{
  return cost >= 0;
}

#include <stdio.h>
#include <ostream>
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <map>
#include <list>
#include "costs.h"
#include <fstream>
#include <cstdlib>
#include "IndStream.h"


using namespace std;

#if DOELEM
enum DistType { UNKNOWN, 
		D_STAR_STAR, 
		D_MC_MR, 
		D_MC_MR_T,
		D_MC_MR_H,
		D_MR_MC, 
		D_MR_MC_T,
		D_MR_MC_H,
		D_MC_STAR, 
		D_MC_STAR_T,
		D_MC_STAR_H,
		D_STAR_MC,
		D_STAR_MC_T,
		D_STAR_MC_H,
		D_MR_STAR, 
		D_MR_STAR_T,
		D_MR_STAR_H,
		D_STAR_MR, 
		D_STAR_MR_T,
		D_STAR_MR_H,
		D_VC_STAR, 
		D_VC_STAR_H,
		D_VC_STAR_T,
		D_STAR_VC, 
		D_VR_STAR, 
		D_STAR_VR,
		D_LASTDIST };
#elif DOTENSORS

#endif

enum Trans { NORMAL,
	     TRANS,
	     CONJTRANS,
	     CONJ };

enum Side { LEFT,
	    RIGHT };

enum Tri { LOWER,
	   UPPER,
	   NOTTRI};

enum Diag { UNIT,
	    NONUNIT,
	    NOTTRIDIAG };

enum TriStruct { HERM,
		 SYMM,
		 TRI,
		 GEN };

enum Type { REAL,
	    COMPLEX };

enum Dir { HORIZONTAL,
	   VERTICAL };

enum PossTunType {
  POSSTUNIN,
  POSSTUNOUT,
  SETTUNIN,
  SETTUNOUT,
  LASTTUNNEL
};

enum PackSize { USEMRSIZE,
		USEKRSIZE,
		USENRSIZE,
		BADPACKSIZE };

enum Dim {
  DIMM,
  DIMK,
  DIMN,
  BADDIM
};

string BoolToStr(bool boolVal);	    
#if DOELEM
string DistTypeToStr(DistType distType);
DistType GetBaseDistType(DistType distType);
#endif
string TransToStr(Trans trans);
char TransToChar(Trans trans);
string SideToStr(Side side);
char SideToChar(Side side);
string TriToStr(Tri tri);
char TriToChar(Tri tri);
Tri SwapTri(Tri tri);
string DiagToStr(Diag diag);

bool IsTrans(Trans trans);

class Node;
class Name;
class Transformation;
class Poss;
class NodeConn;
class PSet;
class InputNode;
typedef string ClassType;
typedef vector<Node*> NodeVec;
typedef NodeVec::iterator NodeVecIter;
typedef NodeVec::const_iterator NodeVecConstIter;
typedef set<Node*> NodeSet;
typedef NodeSet::iterator NodeSetIter;
typedef map<Node*,Node*> NodeMap;
typedef NodeMap::iterator NodeMapIter;
//typedef vector<Symb> SymbVec;
//typedef SymbVec::iterator SymbVecIter;
//typedef SymbVec::const_iterator SymbVecConstIter;
typedef vector<Poss*> PossVec;
typedef PossVec::iterator PossVecIter;
typedef PossVec::const_iterator PossVecConstIter;
typedef vector<PSet*> PSetVec;
typedef PSetVec::iterator PSetVecIter;
typedef PSetVec::const_iterator PSetVecConstIter;
typedef vector<Name> NameVec;
typedef vector<NodeConn*> NodeConnVec;
typedef NodeConnVec::iterator NodeConnVecIter;
typedef NodeConnVec::const_iterator NodeConnVecConstIter;
typedef string NodeType;
typedef vector<Transformation*> TransVec;
typedef TransVec::iterator TransVecIter;
typedef TransVec::const_iterator TransVecConstIter;
typedef vector<const Transformation*> TransConstVec;
typedef TransConstVec::const_iterator TransConstVecIter;
typedef set<const Transformation*> TransSet;
typedef TransSet::iterator TransSetIter;
typedef TransSet::const_iterator TransSetConstIter;
typedef set<const PSet*> PSetSet;
typedef PSetSet::iterator PSetSetIter;
typedef PSetSet::const_iterator PSetSetConstIter;
typedef map<ClassType,TransVec*> TransMap;
typedef TransMap::iterator TransMapIter;
typedef TransMap::const_iterator TransMapConstIter;
typedef vector<string> StrVec;
typedef StrVec::iterator StrVecIter;
typedef StrVec::const_iterator StrVecConstIter;
typedef map<Node*,int> NodeIntMap;
typedef NodeIntMap::iterator NodeIntMapIter;
typedef set<string> StrSet;
typedef StrSet::iterator StrSetIter;
typedef StrSet::const_iterator StrSetConstIter;
typedef set<int> IntSet;
typedef IntSet::iterator IntSetIter;
typedef IntSet::const_iterator IntSetConstIter;
typedef vector<InputNode*> InputVec;
typedef InputVec::iterator InputVecIter;
typedef vector<Transformation*> TransPtrVec;
typedef TransPtrVec::iterator TransPtrVecIter;
typedef TransPtrVec::const_iterator TransPtrVecConstIter;
typedef map<string,const Transformation*> TransNameMap;
typedef TransNameMap::iterator TransNameMapIter;
typedef TransNameMap::const_iterator TransNameMapConstIter;
typedef map<const Transformation*,string> TransPtrMap;
typedef TransPtrMap::iterator TransPtrMapIter;
typedef TransPtrMap::const_iterator TransPtrMapConstIter;
typedef map<void*,const void*> PtrMap;
typedef PtrMap::iterator PtrMapIter;
typedef PtrMap::const_iterator PtrMapConstIter;
//typedef vector<Size> Sizes;
typedef vector<Cost> CostVec;
typedef CostVec::iterator CostVecIter;


bool FoundInNodeVec(const NodeVec &vec, const Node *node);

//Variable name
class Name
{
 public:
#if DOELEM
  DistType m_type;
#endif
  string m_name;
 Name() :
#if DOELEM
  m_type(UNKNOWN), 
#endif
    m_name("noname") {}
  string str() const;
  Name& operator=(const Name &rhs);
  void Flatten(ofstream &out) const;
  void Unflatten(ifstream &in);
};



template<class T>
bool AddElemToVec(std::vector<T*> &vec, T *elem, bool deep = true);
bool AddPossesToVecOrDispose(PossVec &vec, const PossVec &newPoss);


typedef void (*CullFunction)(Poss *poss, bool &cullIfPossible, bool &doNotCull);

struct SaveInfo
{
  PtrMap *transMap;
  PtrMap *possMap;
  PtrMap *psetMap;
  NodeMap *nodeMap;
};

template<class T>
inline void Swap(T **ptr, PtrMap *map)
{
  if (*ptr == NULL)
    return;
  PtrMapIter iter = map->find(*ptr);
  if (iter == map->end()) {
    cout << "didn't find swap pointer\n";
    throw;
  }
  *ptr = (T*)(iter->second);  
}


inline void Swap(Node **ptr, NodeMap *map)
{
  if (*ptr == NULL)
    return;
  NodeMapIter iter = map->find(*ptr);
  if (iter == map->end()) {
    cout << "didn't find swap pointer\n";
    throw;
  }
  *ptr = (Node*)(iter->second);  
}


string LayerNumToStr(Layer layer);

#define SIMP -1
#define GLOBSIMP -2
