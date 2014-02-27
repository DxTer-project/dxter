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



#include "transform.h"
#include "poss.h"
#include <cstring>
#include <cmath>
#include <math.h>

char END = '#';
char START = '+';


template<>
bool AddElemToVec(std::vector<Poss*> &vec, Poss *elem, bool deep)
{
  std::vector<Poss*>::iterator iter = vec.begin();
  for( ; iter != vec.end(); ++iter) {
    if ((deep && **iter == *elem) || (!deep && *iter == elem)) {
      return false;
    }
  }
  vec.push_back(elem);
  return true;
}

template<>
bool AddElemToVec(std::vector<PSet*> &vec, PSet *elem, bool deep)
{
  if (deep)
    throw;
  std::vector<PSet*>::iterator iter = vec.begin();
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

#if DOELEM
string DistTypeToStr(DistType distType)
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
  return m_name + "_" + DistTypeToStr(m_type);
#else
  return m_name;
#endif
}

Name& Name::operator= (const Name &rhs)
{
#if DODM
  m_type = rhs.m_type;
#endif
  m_name = rhs.m_name;
  return *this;
}

void Name::Flatten(ofstream &out) const
{
#if DODM
  WRITE(m_type);
#endif
  out << m_name << endl;
}

void Name::Unflatten(ifstream &in)
{
#if DODM
  READ(m_type);
#endif
  getline(in, m_name);
}

NodeConn::NodeConn(Node *n, unsigned int num)
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


string LayerNumToStr(Layer layer)
{
  switch(layer) {
  case(0):
    return "0";
  case(1):
    return "1";
  case(2):
    return "2";
  case(3):
    return "3";
  case(4):
    return "4";
    case (5):
    return "5";
  default:
    throw;
  }
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
