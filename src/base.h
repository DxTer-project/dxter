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

#include "layers.h"
#include "sizes.h"
#include <cstring>

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

#if TWOD
enum DimName {
  DIMM,
  DIMK,
  DIMN,
  BADDIM
};
#endif

string BoolToStr(bool boolVal);	    

string TransToStr(Trans trans);
char TransToChar(Trans trans);
string SideToStr(Side side);
char SideToChar(Side side);
string TriToStr(Tri tri);
char TriToChar(Tri tri);
Tri SwapTri(Tri tri);
string DiagToStr(Diag diag);

bool IsTrans(Trans trans);


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

extern DistType STAR_STAR; 
extern DistType MC_MR; 
extern DistType MC_MR_T;
extern DistType MC_MR_H;
extern DistType MR_MC; 
extern DistType MR_MC_T;
extern DistType MR_MC_H;
extern DistType MC_STAR; 
extern DistType MC_STAR_T;
extern DistType MC_STAR_H;
extern DistType STAR_MC;
extern DistType STAR_MC_T;
extern DistType STAR_MC_H;
extern DistType MR_STAR; 
extern DistType MR_STAR_T;
extern DistType MR_STAR_H;
extern DistType STAR_MR; 
extern DistType STAR_MR_T;
extern DistType STAR_MR_H;
extern DistType VC_STAR; 
extern DistType VC_STAR_H;
extern DistType VC_STAR_T;
extern DistType STAR_VC; 
extern DistType VR_STAR; 
extern DistType STAR_VR;
#elif DOTENSORS

typedef unsigned int Dim;
typedef vector<Dim> DimVec;
typedef DimVec::iterator DimVecIter;
typedef DimVec::const_iterator DimVecConstIter;
typedef DimVec::reverse_iterator DimVecRevIter;
typedef DimVec::const_reverse_iterator DimVecConstRevIter;
//typedef unordered_set<Dim> DimSet;
typedef set<Dim> DimSet;
typedef DimSet::iterator DimSetIter;
typedef DimSet::const_iterator DimSetConstIter;


class DistEntry
{
 public:
  unsigned int m_val;
  DistEntry() { m_val = 0; }
  inline DistEntry& operator=(const DistEntry &rhs)  { m_val = rhs.m_val; return *this; }
  inline bool operator==(const DistEntry &rhs) const { return m_val == rhs.m_val; }
  inline bool operator!=(const DistEntry &rhs) const { return m_val != rhs.m_val; }
  inline bool operator<(const DistEntry &rhs) const { return m_val < rhs.m_val; }
  inline bool IsStar() const { return m_val == 0; }
  inline void SetToStar() { m_val = 0; }
  string str() const;
  string PrettyStr() const;
  DimVec DistEntryDims() const;
  DimSet DistEntryDimSet() const;
  void DimsToDistEntry(DimVec dims);
};

struct DistEntryCompare {
  bool operator() (const DistEntry& lhs, const DistEntry& rhs) const{
    return lhs < rhs;
  }
};

typedef set<DistEntry,DistEntryCompare> EntrySet;
typedef EntrySet::iterator EntrySetIter;
typedef EntrySet::const_iterator EntrySetConstIter;


class DistType
{
 public:
  Dim m_numDims;
  DistEntry *m_dists;
  DistEntry m_notReped; //not replicated
 DistType() 
   : m_numDims(99), m_dists(NULL) {}
  DistType(const DistType &rhs);
  ~DistType();
  void SetToDefault(Dim numDims);
  void SetToScalarNoRep();
  void SetToStar(Dim numDims);
  DimSet UsedGridDims() const;
  bool IsSane() const;
  void PrepForNumDims(Dim numDims);
  string QuickStr() const;
  string str() const;
  string PrettyStr() const;
  bool HasNoReped() const { return m_notReped.IsStar(); }
  void AddNotReped(Dim dim);
  void AddNotReped(DistEntry entry);
  DistType& operator=(const DistType &rhs);
};
#endif

#if DOELEM
inline bool DistTypeNotEqual(const DistType &one, const DistType &two) 
{
  return one != two;
}
#elif DOTENSORS

inline bool DistTypeNotEqual(const DistType &one, const DistType &two)
{
  if (one.m_numDims != two.m_numDims)
    return true;
  return memcmp(one.m_dists,two.m_dists,one.m_numDims*sizeof(unsigned int )) != 0;
}

inline bool DistTypeEqual(const DistType &one, const DistType &two)
{
  return !DistTypeNotEqual(one, two);
}
#endif

#if DODM
string DistTypeToStr(const DistType &type);
void CheckThisSpot();
#if TWOD
void GetLocalSizes(DistType dist, const Sizes *m, const Sizes *n, Sizes &localM, Sizes &localN);
#else
void GetLocalSizes(DistType dist, const SizesArray sizes, SizesArray localSizes);
void GetLocalSizes(DistType dist, Dim dim, const Sizes* sizes, Sizes* localSizes);
#endif
#endif

#if DOELEM
DistType GetBaseDistType(DistType distType);
#endif


class Node;
class Name;
class Transformation;
class Poss;
class NodeConn;
//class BasePSet;
class InputNode;
typedef unsigned long int GraphNum;
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
typedef multimap<size_t,Poss*> PossMMap;
typedef PossMMap::iterator PossMMapIter;
typedef PossMMap::const_iterator PossMMapConstIter;
typedef std::pair<size_t,Poss*> PossMMapPair;
typedef std::pair<PossMMapIter,PossMMapIter> PossMMapRangePair;
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
typedef vector<double> TimeVec;
typedef TimeVec::iterator TimeVecIter;
typedef TimeVec::const_iterator TimeVecConstIter;
typedef map<GraphNum, TimeVec> ImplementationRuntimeMap;
typedef ImplementationRuntimeMap::iterator ImplementationRuntimeMapIter;
typedef ImplementationRuntimeMap::const_iterator ImplementationRuntimeMapConstIter;

double MinTime(const TimeVec &times);



#if DOTENSORS
typedef vector<DistType*> DistTypeVec;
typedef DistTypeVec::iterator DistTypeVecIter;
typedef DistTypeVec::const_iterator DistTypeVecConstIter;
DimVec MapIndicesToDims(const string &indices, const string &dimIndices);
bool IsPrefix(const DimVec &isPrefix, const DimVec &dims);
#endif



//Variable name
class Name
{
 public:
#if DOELEM
  DistType m_type;
#elif DOTENSORS
  DistType m_type;
#endif
  string m_name;
 Name() :
#if DOELEM
  m_type(UNKNOWN), 
#endif
    m_name("noname") {}
  string str() const;
#if DOTENSORS
  string PrettyStr() const;
#endif
  Name& operator=(const Name &rhs);
  void Flatten(ofstream &out) const;
  void Unflatten(ifstream &in);
};



bool FoundInNodeVec(const NodeVec &vec, const Node *node);



template<class T>
bool AddElemToVec(std::vector<T*> &vec, T *elem, bool deep = true);

bool AddPossToMMap(PossMMap &mmap, Poss *elem, size_t hash, bool deep = true);


//bool AddPossesToVecOrDispose(PossVec &vec, const PossVec &newPoss);


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


inline string LayerNumToStr(Layer layer) 
{ 
  return std::to_string(layer); 
}

#define SIMP -1
#define GLOBSIMP -2



