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



#include "base.h"
#include <cmath>
#include "transform.h"
#include "poss.h"
#include <cstring>
#include <math.h>
#include "costs.h"
#include <ostream>
#include <sstream>

using namespace std;

#if DOELEM
DistType STAR_STAR = D_STAR_STAR;
DistType MC_MR = D_MC_MR;
DistType MC_MR_T = D_MC_MR_T;
DistType MC_MR_H = D_MC_MR_H;
DistType MR_MC = D_MR_MC;
DistType MR_MC_T = D_MR_MC_T;
DistType MR_MC_H = D_MR_MC_H;
DistType MC_STAR = D_MC_STAR;
DistType MC_STAR_T = D_MC_STAR_T;
DistType MC_STAR_H = D_MC_STAR_H;
DistType STAR_MC = D_STAR_MC;
DistType STAR_MC_T = D_STAR_MC_T;
DistType STAR_MC_H = D_STAR_MC_H;
DistType MR_STAR = D_MR_STAR;
DistType MR_STAR_T = D_MR_STAR_T;
DistType MR_STAR_H = D_MR_STAR_H;
DistType STAR_MR = D_STAR_MR;
DistType STAR_MR_T = D_STAR_MR_T;
DistType STAR_MR_H = D_STAR_MR_H;
DistType VC_STAR = D_VC_STAR;
DistType VC_STAR_H = D_VC_STAR_H;
DistType VC_STAR_T = D_VC_STAR_T;
DistType STAR_VC = D_STAR_VC;
DistType VR_STAR = D_VR_STAR;
DistType STAR_VR = D_STAR_VR;
#endif

char END = '#';
char START = '+';


#if DOTENSORS
DistType::DistType(const DistType &rhs)
  :m_dists(NULL)
{
  m_numDims = rhs.m_numDims;
  if (m_numDims) {
    m_dists = new DistEntry[m_numDims];
    for (Dim dim = 0; dim < m_numDims; ++dim)
      m_dists[dim] = rhs.m_dists[dim];
  }
  m_notReped = rhs.m_notReped;
}

DistType& DistType::operator=(const DistType &rhs)
{
  m_numDims = rhs.m_numDims;
  m_notReped = rhs.m_notReped;
  if (m_dists)
    delete [] m_dists;
 if (m_numDims) {
    m_dists = new DistEntry[m_numDims];
    for (Dim dim = 0; dim < m_numDims; ++dim)
      m_dists[dim] = rhs.m_dists[dim];
  }
  else
    m_dists = NULL;
  return *this;
}

DistType DistType::Permute(const Permutation &perm) const
{
  DistType ret;
  ret.m_notReped = m_notReped;
  ret.PrepForNumDims(m_numDims);
  for(Dim dim = 0; dim < m_numDims; ++dim) {
    ret.m_dists[perm.MapStartToFinish(dim)] = m_dists[dim];
  }
  return ret;
}

void DistType::PrepForNumDims(Dim numDims)
{
  if (m_dists)
    delete [] m_dists;
  m_numDims = numDims;
  if (numDims)
    m_dists = new DistEntry[numDims];
  else
    m_dists = NULL;
  m_notReped.SetToStar();
}

DistType::~DistType()
{
  if (m_dists) {
    delete [] m_dists;
    m_dists = NULL;
  }
  m_numDims = 0;
}

DimSet DistType::UsedGridDims() const
{
  //Reflect changes in IsSane
  DimSet set;
  for (Dim dim = 0; dim < m_numDims; ++dim) {
    DimVec vec = m_dists[dim].DistEntryDims();
    DimVecIter iter = vec.begin();
    for(; iter != vec.end(); ++iter) {
      if (!set.insert(*iter).second) {
	cout << DistTypeToStr(*this) << endl;
	throw;
      }
    }
  }
  DimSet repedSet = m_notReped.DistEntryDimSet();
  set.insert(repedSet.begin(), repedSet.end());
  return set;
}

bool DistType::IsSane() const
{
  //Reflect changes in UsedGridDims
  DimSet set;
  for (Dim dim = 0; dim < m_numDims; ++dim) {
    DimVec vec = m_dists[dim].DistEntryDims();
    DimVecIter iter = vec.begin();
    for(; iter != vec.end(); ++iter) {
      if (!set.insert(*iter).second) {
	return false;
      }
    }
  }
  DimVec vec = m_notReped.DistEntryDims();
  DimVecIter iter = vec.begin();
  for(; iter != vec.end(); ++iter) {
    if (!set.insert(*iter).second) {
      return false;
      }
  }
  return true;
}


void DistType::SetToDefault(Dim numDims)
{
  if (numDims > NUM_GRID_DIMS)
    throw;
  m_numDims = numDims;
  if (m_dists)
    delete [] m_dists;
  m_notReped.SetToStar();
  m_dists = new DistEntry[numDims];

  unsigned int numStartDists = ceil((double)NUM_GRID_DIMS / numDims);
  unsigned int numEndDists = floor((double)NUM_GRID_DIMS / numDims);
  unsigned int tipPoint = numDims;
  if (numEndDists != numStartDists) 
    tipPoint = (NUM_GRID_DIMS - numEndDists * numDims) / (numStartDists - numEndDists);

  /*    
  cout <<"numDims: " << numDims 
       << "\ngridDims: " << NUM_GRID_DIMS << endl;
  
  cout << "start " << numStartDists << endl;
  cout << "end " << numEndDists << endl;
  cout << "tip " << tipPoint << endl;
  */
  Dim currDistDim = 0;
  for (Dim i = 0; i < numDims; ++i) {
    DimVec vec;
    for (unsigned int j = 0; j < (i < tipPoint ? numStartDists : numEndDists); ++j)  {
      vec.push_back(currDistDim);
      currDistDim++;
    }
    m_dists[i].DimsToDistEntry(vec);
  }
}

void DistType::SetToScalarNoRep()
{
  m_numDims = 0;
  if (m_dists) {
    delete [] m_dists;
    m_dists = NULL;
  }
  m_notReped.SetToStar();
  DimVec vec;
  for (Dim j = 0; j < NUM_GRID_DIMS; ++j)  {
      vec.push_back(j);
  }
    m_notReped.DimsToDistEntry(vec);
}

void DistType::SetToStar(Dim numDims)
{
  if (numDims > NUM_GRID_DIMS)
    throw;
  m_numDims = numDims;
  if (m_dists)
    delete [] m_dists;
  m_dists = new DistEntry[numDims];
  m_notReped.SetToStar();

  for(Dim dim = 0; dim < numDims; ++dim)
    m_dists[dim].SetToStar();
}

string DistEntry::str() const
{
    DimVec vec = DistEntryDims();
  if (vec.empty())
    return "S";
  std::stringstream ret;
  ret << "D";
  DimVecIter iter = vec.begin();
  for (; iter != vec.end(); ++iter) {
    ret << "_";
    ret << *iter;
  }
  return ret.str();
}


string DistEntry::PrettyStr() const
{
  DimVec vec = DistEntryDims();
  if (vec.empty())
    return "*";
  std::stringstream ret;
  ret << "D";
  DimVecIter iter = vec.begin();
  for (; iter != vec.end(); ++iter) {
    ret << *iter;
  }
  return ret.str();
}



string DistType::QuickStr() const
{
  if (!m_numDims) 
    return "";
  std::stringstream ret;
  if (!m_dists)
    throw;
  for(Dim dim = 0; dim < m_numDims; ++dim) {
    ret << m_dists[dim].m_val << " ";
  }
  ret << "|" << m_notReped.m_val;    
  return ret.str();
}

DimVec DistEntry::DistEntryDims() const
{
#if 0
  if (dist > 3000) {
    cout << "big dist to start: " << dist << endl;
    throw;
  }
  cout << "\n\nstarting with dist = " << dist << endl;
#endif
  DimVec vec;
  if (IsStar())
    return vec;
  unsigned int currStage = NUM_GRID_DIMS;
  unsigned int distVal = m_val-1;
  unsigned int numDists = 1;
  while (distVal > currStage) {
  //  while (m_val > currStage) {
    if (!currStage) {
      cout << "starting dist " << m_val << endl;
      cout << "numDists " << numDists << endl;
      cout << "!currStage\n";
      throw;
    }
#if 0
    cout << "distVal bef " << distVal << endl;
    cout << "numDists bef " << numDists << endl;
    cout << "currStage bef " << currStage << endl;
    if ((distVal - currStage) >= distVal) {
      cout << "problem with current\n";
      throw;
    }
#endif
    distVal -= currStage;
    numDists++;
    currStage *= NUM_GRID_DIMS;
#if 0
    cout << "distVal aft " << distVal << endl;
    cout << "numDists aft " << numDists << endl;
    cout << "currStage aft " << currStage << endl;
#endif
  }
  string out;
  while (numDists) {
#if 0
    cout << "inserting " << distVal % MAX_NUM_DIMS << endl;
    cout << "overwriting " << distVal << " with " << distVal / MAX_NUM_DIMS << endl;
    cout << "numDists " << numDists << " being decremented" << endl;
#endif
    vec.insert(vec.begin(), distVal % NUM_GRID_DIMS);
    distVal = distVal / NUM_GRID_DIMS;
    --numDists;
  }
  if (distVal != 0) {
    cout << endl << distVal << " != 0\n";
    cout << "dist = " << m_val << endl;
    throw;
  }
  return vec;
}

DimSet DistEntry::DistEntryDimSet() const
{
  DimVec vec(DistEntryDims());
  DimSet set;
  set.insert(vec.begin(), vec.end());
  return set;
}

void DistEntry::DimsToDistEntry(const DimVec &dims)
{
  unsigned int currStage = 1;
  unsigned int distVal = 0;
  DimVecConstRevIter iter = dims.rbegin();
  for(; iter != dims.rend(); ++iter) {
    if (*iter > NUM_GRID_DIMS) {
      cout << *iter << endl;
      throw;
    }
    //offset coming from the left
    distVal += currStage;
    //offset on the right
    distVal += currStage * *iter;
    
    currStage *= NUM_GRID_DIMS;
  }

  if (distVal > 1410273309) {
    cout << "big distVal: " << distVal << endl;
    throw;
  }
    
  m_val = distVal;
}

string DistType::str() const
{
  string out = "_";
  for (Dim i = 0; i < m_numDims; ++i) {
    out += m_dists[i].str();
    if (i+1 < m_numDims)
      out += "__";
  }
  if (!m_notReped.IsStar()) {
    out += "__N_";
    out += m_notReped.str();
  }
  return out;
}

string DistType::PrettyStr() const
{
  std::stringstream out;
  out << "[";
  for (Dim i = 0; i < m_numDims; ++i) {
    out << m_dists[i].PrettyStr();
    if (i+1 < m_numDims)
      out << ",";
  }
  out << "]";
  if (!m_notReped.IsStar()) {
    out << " | {";
    DimVec vec = m_notReped.DistEntryDims();
    DimVecIter iter = vec.begin();
    out << *iter;
    ++iter;
    for(; iter != vec.end(); ++iter) {
      out << "," << *iter;
    }
    out << "}";      
  }

  return out.str();
}

void DistType::AddNotReped(Dim dim)
{
  DimSet set = m_notReped.DistEntryDimSet();
  set.insert(dim);
  DimVec vec;
  vec.insert(vec.end(), set.begin(), set.end());
  m_notReped.DimsToDistEntry(vec);
}

void DistType::AddNotReped(DistEntry entry)
{
  DimSet set = m_notReped.DistEntryDimSet();
  DimSet set2 = entry.DistEntryDimSet();
  set.insert(set2.begin(),set2.end());
  DimVec vec;
  vec.insert(vec.end(), set.begin(), set.end());
  m_notReped.DimsToDistEntry(vec);
}



string DistTypeToStr(const DistType &type)
{
  return type.str();
}

#if DOTENSORS
DimVec MapIndicesToDims(const string &indices, const string &dimIndices)
{
  DimVec map;
  map.reserve(indices.length());
  string::const_iterator iter = indices.begin();
  for(; iter != indices.end(); ++iter) {
    map.push_back(dimIndices.find(*iter));
  }
  return map;
}

bool IsPrefix(const DimVec &isPrefix, const DimVec &dims)
{
  if (dims.empty() || isPrefix.empty())
    return false;
  DimVecConstIter dimsIter = dims.begin();
  DimVecConstIter isPrefixIter = isPrefix.begin();
  for(; dimsIter != dims.end() && isPrefixIter != isPrefix.end(); ++dimsIter, ++isPrefixIter) {
    if (*isPrefixIter != *dimsIter)
      return false;
  }
  if (isPrefixIter != isPrefix.end() && dimsIter == dims.end())
    return false;
  return true;
}
#endif

void GetLocalSizes(const DistType &dist, const SizesArray sizes, SizesArray localSizes)
{
  const Dim numDims = dist.m_numDims;
  for (Dim dim = 0; dim < numDims; ++ dim) {
    const DistEntry &entry = dist.m_dists[dim];
    localSizes[dim] = sizes[dim];
    if (!entry.IsStar()) {
      DimVec vec = entry.DistEntryDims();
      unsigned int coef = 1;
      DimVecIter iter = vec.begin();
      for(; iter != vec.end(); ++iter) {
	coef *= GridLens[*iter];
      }
      localSizes[dim].SetCoeff(1.0 / coef);
    }
  }
}

void GetLocalSizes(const DistType &dist, Dim dim, const Sizes* sizes, Sizes* localSizes)
{
  const DistEntry &entry = dist.m_dists[dim];
  *localSizes = *sizes;
  if (!entry.IsStar()) {
    DimVec vec = entry.DistEntryDims();
    unsigned int coef = 1;
    DimVecIter iter = vec.begin();
    for(; iter != vec.end(); ++iter) {
      coef *= GridLens[*iter];
    }
    localSizes->SetCoeff(1.0 / coef);
  }
}
#endif

#if DOELEM
string DistTypeToStr(const DistType &distType)
{
  switch(distType) {
    case D_STAR_STAR:
      return "STAR_STAR";
      break;
    case D_MC_MR:
      return "MC_MR";
      break;
    case D_VC_STAR_H:
      return "Adj_VC_STAR";
      break;
    case D_VC_STAR_T:
      return "Trans_VC_STAR";
      break;
    case D_MC_MR_T:
      return "Trans_MC_MR";
      break;
    case D_MC_MR_H:
      return "Adj_MC_MR";
      break;
    case D_MR_MC:
      return "MR_MC";
      break;
    case D_MR_MC_T:
      return "Trans_MR_MC";
      break;
    case D_MR_MC_H:
      return "Adj_MR_MC";
      break;
    case D_MC_STAR:
      return "MC_STAR";
      break;
    case D_MC_STAR_T:
      return "Trans_MC_STAR";
      break;
    case D_MC_STAR_H:
      return "Adj_MC_STAR";
      break;
    case D_STAR_MC:
      return "STAR_MC";
      break;
    case D_STAR_MC_T:
      return "Trans_STAR_MC";
      break;
    case D_STAR_MC_H:
      return "Adj_STAR_MC";
      break;
    case D_MR_STAR:
      return "MR_STAR";
      break;
    case D_MR_STAR_T:
      return "Trans_MR_STAR";
      break;
    case D_MR_STAR_H:
      return "Adj_MR_STAR";
      break;
    case D_STAR_MR:
      return "STAR_MR";
      break;
    case D_STAR_MR_T:
      return "Trans_STAR_MR";
      break;
    case D_STAR_MR_H:
      return "Adj_STAR_MR";
      break;
    case D_VC_STAR:
      return "VC_STAR";
      break;
    case D_STAR_VC:
      return "STAR_VC";
      break;
    case D_VR_STAR:
      return "VR_STAR";
      break;
    case D_STAR_VR:
      return "STAR_VR";
      break;
    default:
      throw;
      break;
  }
}


DistType GetBaseDistType(DistType distType)
{
  switch(distType) {
    case D_MC_STAR_T:
    case D_MC_STAR_H:
      return D_MC_STAR;
    case D_STAR_MC_T:
    case D_STAR_MC_H:
      return D_STAR_MC;
    case D_MR_STAR_T:
    case D_MR_STAR_H:
      return D_MR_STAR;
    case D_STAR_MR_T:
    case D_STAR_MR_H:
      return D_STAR_MR;
    default:
      return distType;
  }
}


#endif //DOELEM



bool AddPossToMMap(PossMMap &mmap, Poss *elem, size_t hash, bool deep)
{
  PossMMapRangePair pair = mmap.equal_range(hash);
  for( ; pair.first != pair.second; ++pair.first) {
    Poss *poss = (*(pair.first)).second;
    if ((deep && *poss == *elem) || (!deep && poss == elem)) {
      return false;
    }
  }
  mmap.insert(PossMMapPair(hash, elem));
  return true;
}

template<>
bool AddElemToVec(std::vector<BasePSet*> &vec, BasePSet *elem, bool deep)
{
  if (deep)
    throw;
  std::vector<BasePSet*>::iterator iter = vec.begin();
  for( ; iter != vec.end(); ++iter) {
    if (*iter == elem) {
      return false;
    }
  }
  vec.push_back(elem);
  return true;
}

template<class T>
bool AddElemToVec(std::vector<T*> &vec, T *elem, bool deep)
{
  typename std::vector<T*>::iterator iter = vec.begin();
  for( ; iter != vec.end(); ++iter) {
    if ((deep && **iter == *elem) || (!deep && *iter == elem)) {
      return false;
    }
  }
  vec.push_back(elem);
  return true;
}

template bool AddElemToVec<Node>(std::vector<Node*> &vec, Node *elem, bool deep);
template bool AddElemToVec<Poss>(std::vector<Poss*> &vec, Poss *elem, bool deep);


string BoolToStr(bool boolVal)
{
  if (boolVal)
    return "TRUE";
  else
    return "FALSE";
}

string DiagToStr(Diag diag)
{
  if (diag == UNIT)
    return "UNIT";
  else if (diag == NONUNIT)
    return "NONUNIT";
  else
    throw;
}


string TransToStr(Trans trans)
{
  switch(trans) {
    case (NORMAL):
      return "NORMAL";
    case (TRANS):
      return "TRANSPOSE";
    case (CONJTRANS):
      return "ADJOINT";
    default: {
      return "unknown";
    }
  }
}

char TransToChar(Trans trans)
{
  switch(trans) {
    case (NORMAL):
      return 'N';
    case (TRANS):
      return 'T';
    case (CONJTRANS):
      return 'H';
    default:
      return 'U';
  }
}


string SideToStr(Side side)
{
  switch(side) {
    case (LEFT):
      return "LEFT";
    case (RIGHT):
      return "RIGHT";
    default:
      return "unknown";
  }
}

char SideToChar(Side side)
{
  switch(side) {
    case (LEFT):
      return 'L';
    case (RIGHT):
      return 'R';
    default:
      return 'U';
  }
}

string TriToStr(Tri tri)
{
  switch(tri) {
    case (LOWER):
      return "LOWER";
    case (UPPER):
      return "UPPER";
    default:
      return "unknown";
  }
}

char TriToChar(Tri tri)
{
  switch(tri) {
    case (LOWER):
      return 'L';
    case (UPPER):
      return 'U';
    default:
      return 'K';
  }
}

Tri SwapTri(Tri tri)
{
  switch(tri) {
  case(LOWER):
    return UPPER;
  case(UPPER):
    return LOWER;
  default:
    throw;
  }
}

bool IsTrans(Trans trans)
{
  return trans == TRANS || trans == CONJTRANS;
}


string Name::str() const
{
#if DOELEM
  string name = m_name;
  if (m_type == D_MC_MR)
    return name;
  else if (m_type == UNKNOWN)
    return name;
  else
    return name + "_" + DistTypeToStr(m_type);
#elif DOTENSORS
  string name = m_name;
  if (m_permutation.Size()) {
    name += "_perm" + m_permutation.Str();
    name += "_" + DistTypeToStr(m_type.Permute(m_permutation));
  }
  else {
    name += "_" + DistTypeToStr(m_type);
  }
  return name;
#else
  return m_name;
#endif
}

#if DOTENSORS
string Name::PrettyStr() const
{
  return m_name + m_type.PrettyStr();
}
#endif

Name& Name::operator= (const Name &rhs)
{
#if DODM
  m_type = rhs.m_type;
#endif
#if DOTENSORS
  m_permutation = rhs.m_permutation;
#endif
  m_name = rhs.m_name;
  return *this;
}

void Name::Flatten(ofstream &out) const
{
#if DODM
  WRITE(m_type);
#endif
#if DOTENSORS
  throw;
  //m_permutation
#endif
  out << m_name << endl;

}

void Name::Unflatten(ifstream &in)
{
#if DODM
  READ(m_type);
#endif
#if DOTENSORS
  throw;
  //m_permutation
#endif
  getline(in, m_name);
}

NodeConn::NodeConn(Node *n, ConnNum num)
: m_n(n), m_num(num)
{
  if (!n)
    cout << "!n\n";
}

NodeConn::NodeConn(const NodeConn *conn)
: m_n(conn->m_n), m_num(conn->m_num)
{
  if (!m_n)
    cout << "!m_n\n";
}

bool NodeConn::operator==(const NodeConn &rhs) const
{
  if (m_num != rhs.m_num)
    return false;
  if (!m_n || !rhs.m_n) {
    cout << "bad conn\n";
    throw;
  }
  return *m_n == *(rhs.m_n);
}

void NodeConn::SetNode(Node *node)
{
  if (!node) {
    cout << "SetNode(null)\n";
    throw;
  }
  m_n = node;
}

void NodeConn::Flatten(ofstream &out) const
{
  WRITE(m_n);
  WRITE(m_num);
}

void NodeConn::Unflatten(ifstream &in)
{
  READ(m_n);
  READ(m_num);
}

bool FoundInNodeVec(const NodeVec &vec, const Node *node)
{
  NodeVecConstIter iter = vec.begin();
  for(; iter != vec.end(); ++iter)
    if (*iter == node)
      return true;
  return false;
}

unsigned int FindInNodeVec(const NodeVec &vec, const Node *node)
{
  unsigned int i = 0;
  NodeVecConstIter iter = vec.begin();
  for(; iter != vec.end(); ++iter,++i)
    if (*iter == node)
      return i;
  throw;
}

unsigned int FindInSetVec(const PSetVec &vec, const BasePSet *set)
{
  unsigned int i = 0;
  PSetVecConstIter iter = vec.begin();
  for(; iter != vec.end(); ++iter,++i)
    if (*iter == set)
      return i;
  throw;
}


#if DOELEM
void GetLocalSizes(DistType dist, const Sizes *m, const Sizes *n, Sizes &localM, Sizes &localN)
{
  switch(dist) {
  case D_STAR_STAR:
    localM = *m;
    localN = *n;
    break;
  case D_MC_MR:
    localM = *m;
    localM.SetCoeff(1.0 / R);
    localN = *n;
    localN.SetCoeff(1.0 / C);
    break;
  case D_MR_MC:
    localM = *m;
    localM.SetCoeff(1.0/C);
    localN = *n;
    localN.SetCoeff(1.0/R);
      break;
    case D_MC_STAR:
      localM = *m;
      localM.SetCoeff(1.0/R);
      localN = *n;
      break;
    case D_MC_STAR_T:
    case D_MC_STAR_H:
      localM = *n;
      localM.SetCoeff(1.0/R);
      localN = *m;
      break;
    case D_STAR_MC:
      localM = *m;
      localN = *n;
      localN.SetCoeff(1.0/R);
      break;
    case D_STAR_MC_T:
    case D_STAR_MC_H:
      localM = *n;
      localN = *m;
      localN.SetCoeff(1.0/ R);
      break;
    case D_MR_STAR:
      localM = *m;
      localM.SetCoeff(1.0 / C);
      localN = *n;
      break;
    case D_MR_STAR_T:
    case D_MR_STAR_H:
      localM = *n;
      localM.SetCoeff(1.0 / C);
      localN = *m;
      break;
    case D_STAR_MR:
      localM = *m;
      localN = *n;
      localN.SetCoeff(1.0 / C);
      break;
    case D_STAR_MR_T:
    case D_STAR_MR_H:
      localM = *n;
      localN = *m;
      localN.SetCoeff(1.0 / C);
      break;
    case D_VC_STAR:
    case D_VR_STAR:
      localM = *m;
      localM.SetCoeff(1.0 / P);
      localN = *n;
      break;
    case D_STAR_VC:
    case D_STAR_VR:
      localM = *m;
      localN = *n;
      localN.SetCoeff(1.0 / P);
      break;
    case D_MC_MR_H:
    case D_MC_MR_T:
      localM = *n;
      localM.SetCoeff(1.0 / R);
      localN = *m;
      localN.SetCoeff(1.0 / C);
      break;
    case D_MR_MC_T:
    case D_MR_MC_H:
      localM = *n;
      localM.SetCoeff(1.0 / C);
      localN = *m;
      localN.SetCoeff(1.0 / R);
    case D_VC_STAR_T:
    case D_VC_STAR_H:
      localM = *n;
      localM.SetCoeff(1.0 / P);
      localN = *m;
      break;
    default:
      throw;
  }
}
#endif


double MinTime(const TimeVec &times)
{
  double minVal = -1;
  TimeVecConstIter iter = times.begin();
  for(; iter != times.end(); ++iter) {
    if (minVal < 0)
      minVal = *iter;
    else if (minVal > *iter)
      minVal = *iter;
  }
  return minVal;
}
