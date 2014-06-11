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



#include "pset.h"

PossTunnel::PossTunnel() 
 :m_tunType(LASTTUNNEL),
  m_pset(NULL)
{
}

PossTunnel::PossTunnel(PossTunType type) 
 :m_tunType(type),
  m_pset(NULL)
{
}

NodeType PossTunnel::GetType() const 
{
  switch(m_tunType)
    {
    case(POSSTUNIN):
      return "Poss input";
    case(POSSTUNOUT):
      return "Poss output";
    case(SETTUNIN):
      return "PSet input";
    case(SETTUNOUT):
      return "PSet output";
    default:
      return "PossTunnel type unknown";
    }
}

bool PossTunnel::IsPossTunnel(PossTunType type) const
{
  return m_tunType == type;
}



void PossTunnel::SetPSet(PSet *set)
{
  if (m_tunType != SETTUNIN && m_tunType != SETTUNOUT)
    cout << "bad set\n";
  m_pset = set;
}

bool PossTunnel::KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const
{
  if (Input(0) == input && InputConnNum(0) == numIn) {
    numOut = numIn;
    return true;
  }
  else
    throw;
}

#if DODM
const DistType& PossTunnel::GetDistType(unsigned int num) const 
{

  if (m_tunType == SETTUNIN || m_tunType == POSSTUNOUT) {
    if (num >= m_inputs.size()) {
      if (m_tunType == SETTUNIN)
	cout << "SETTUNIN\n";
      else
	cout << "POSSTUNOUT\n";
      cout << "num >= m_inputs.size()\n";
      cout << "num = " << num << " and m_inputs.size() = " << m_inputs.size() << endl;
      throw;
    }
    return InputDistType(num); 
  }
  else if (m_tunType == POSSTUNIN || m_tunType == SETTUNOUT) {
    return ((DLANode*)Input(0))->GetDistType(num);
  }
  else {
    cout << "bad tun type\n";
    throw;
  }
}
#endif

PossTunnel* PossTunnel::GetSetTunnel()
{
  if (m_tunType == POSSTUNIN) {
    PossTunnel *tun = (PossTunnel*)GetNewInst();
    tun->m_tunType = SETTUNIN;
    return tun;
  }
  else if (m_tunType == POSSTUNOUT) {
    PossTunnel *tun = (PossTunnel*)GetNewInst();
    tun->m_tunType = SETTUNOUT;
    return tun;
  }
  else {
    cout << "GetSetTunnel on wrong tunnel type\n";
    throw;
  }
}

void PossTunnel::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLANode::Prop();

    if ((m_tunType == SETTUNIN || m_tunType == SETTUNOUT) && !m_pset)
      throw;
    if ((m_tunType == POSSTUNIN && m_inputs.size() != 1)
	|| (m_tunType == POSSTUNOUT && m_children.size() != 1))
      {
	cout << "m_inputs.size() != 1\n";
	if (m_tunType == POSSTUNIN)
	  cout << "m_tunType == POSSTUNIN\n";
	else
	  cout << "m_tunType == POSSTUNOUT\n";
	throw;
      }
    if (m_tunType == LASTTUNNEL) {
      cout << "bad tunnel type!\n";
      throw;
    }
    if (m_tunType == SETTUNIN || m_tunType == POSSTUNIN)
      if (!m_inputs.size()) {
	cout << "!m_inputs.size() on this = " << this << "\n";
	if (m_tunType == POSSTUNIN)
	  cout << "m_tunType == POSSTUNIN\n";
	else 
	  cout << "m_tunType == SETTUNIN is on " << m_pset << endl;
	throw;
      }


    if (m_tunType == POSSTUNIN) {
      if (m_inputs.size() != 1) {
	cout << "m_inputs.size() != 1\n";
	throw;
      }
    }    
    else if (m_tunType == POSSTUNOUT) {
      if (!IsLoopTunnel() && m_inputs.size() != 1) {
	throw;
      }
    }
    if (m_tunType == SETTUNIN || m_tunType == POSSTUNIN)
      if (!m_inputs.size()) {
	cout << "!m_inputs.size() on " << this << "\n";
	throw;
      }

    m_cost = ZERO;
    if (m_tunType == SETTUNOUT) {
      m_pset->Prop();
    }
    else if (m_tunType == POSSTUNOUT || m_tunType == SETTUNIN) {
      for (unsigned int i = 0; i < m_inputs.size(); ++i) {
	Input(i)->Prop();
      }
    }
  }
}

#if TWOD
const Sizes* PossTunnel::GetM(unsigned int num) const
{
  if (m_tunType != SETTUNOUT && m_tunType != POSSTUNIN) {
    return GetInputM(num);
  }
  else {
    return ((DLANode*)Input(0))->GetM(num);
  }
}

const Sizes* PossTunnel::GetN(unsigned int num) const
{
  if (m_tunType != SETTUNOUT && m_tunType != POSSTUNIN) {
    return GetInputN(num);
  }
  else {
    return ((DLANode*)Input(0))->GetN(num);
  }
}

const Sizes* PossTunnel::LocalM(unsigned int num) const
{
  const NodeConn *conn = m_inputs[0];
  const DLANode *input = (DLANode*)(conn->m_n);
  
  return input->LocalM(conn->m_num);
}

const Sizes* PossTunnel::LocalN(unsigned int num) const
{
  const NodeConn *conn = m_inputs[0];
  const DLANode *input = (DLANode*)(conn->m_n);
  
  return input->LocalN(conn->m_num);
}

#if DOLLDLA
Stride PossTunnel::RowStride(unsigned int num) const
{
  return InputRowStride(0);
}

Stride PossTunnel::ColStride(unsigned int num) const
{
  return InputColStride(0);
}
#endif //DOLLDLA

#else
const Dim PossTunnel::NumDims(unsigned int num) const
{
  if (m_tunType != SETTUNOUT && m_tunType != POSSTUNIN) {
    return InputNumDims(num);
  }
  else {
    return ((DLANode*)Input(0))->NumDims(num);
  }
}

const Sizes* PossTunnel::Len(unsigned int num,Dim dim) const
{
  if (m_tunType != SETTUNOUT && m_tunType != POSSTUNIN) {
    return InputLen(num,dim);
  }
  else {
    return ((DLANode*)Input(0))->Len(num,dim);
  }
}

const Sizes* PossTunnel::LocalLen(unsigned int num,Dim dim) const
{
  const NodeConn *conn = m_inputs[0];
  const DLANode *input = (DLANode*)(conn->m_n);
  
  return input->LocalLen(conn->m_num,dim);
}

#endif

Name PossTunnel::GetName(unsigned int num) const 
{
  if (m_tunType != SETTUNOUT && m_tunType != POSSTUNIN) {
    return GetInputName(num);
  }
  else {
    if (num > 0)
      throw;
    return GetInputName(0);
  }
}

unsigned int PossTunnel::NumOutputs() const
{
  switch(m_tunType) 
    {
    case (POSSTUNOUT) :
    case (SETTUNIN) :
      return m_inputs.size();
    case (POSSTUNIN) :
    case (SETTUNOUT) :
	return Input(0)->NumOutputs();
    default:
      cout << "bad tunnel type\n";
      throw;
      return -1;
  } 
}

void PossTunnel::Duplicate(const Node *node, bool shallow, bool possMerging)
{
  DLANode::Duplicate(node, shallow, possMerging);
  m_tunType = ((PossTunnel*)node)->m_tunType;
}

string TunTypeToStr(PossTunType type)
{
  switch(type)
    {
    case (POSSTUNIN):
      return "poss tun in";
    case (POSSTUNOUT):
      return "poss tun out";
    case (SETTUNIN):
      return "set tun in";
    case (SETTUNOUT):
      return "set tun out";
    default:
      return "unknown tun type";
      
    }
}



void PossTunnel::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  WRITE(m_tunType);
  WRITE(m_pset);
}



void PossTunnel::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  DLANode::UnflattenCore(in, info);
  READ(m_tunType);
  READ(m_pset);
  Swap(&m_pset,info.psetMap);
}

bool PossTunnel::Overwrites(const Node *input, unsigned int num) const
{
  if (m_tunType == SETTUNIN) {
    if (!m_children.size())
      throw;
    const NodeConn *conn = m_children[0];
    return conn->m_n->Overwrites(this, conn->m_num);
  }
  else if (m_tunType == POSSTUNIN) {
    NodeConnVecConstIter iter = m_children.begin();
    for(; iter != m_children.end(); ++iter) {
      const NodeConn *conn = *iter;
      if (conn->m_n->Overwrites(this, conn->m_num))
	return true;
    }
    return false;
  }
  else
    return false;
}
