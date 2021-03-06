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


#include "tunnel.h"
#include "basePSet.h"
#include "realPSet.h"
#include "costs.h"

Tunnel::Tunnel() 
 :m_tunType(LASTTUNNEL),
  m_pset(NULL)
#if !DODLA
 ,m_cost(-1)
#endif
{
}

Tunnel::Tunnel(TunType type) 
 :m_tunType(type),
  m_pset(NULL)
#if !DODLA
 ,m_cost(-1)
#endif
{
}

NodeType Tunnel::GetType() const 
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
      return "Tunnel type unknown";
    }
}

bool Tunnel::IsTunnel(TunType type) const
{
  return m_tunType == type;
}


void Tunnel::SetPSet(BasePSet *set)
{
  if (m_tunType != SETTUNIN && m_tunType != SETTUNOUT)
    cout << "bad set\n";
  m_pset = set;
}
/*
bool Tunnel::KeepsInputVarLive(Node *input, ConnNum numIn, ConnNum &numOut) const
{
  if (Input(0) == input && InputConnNum(0) == numIn) {
    numOut = numIn;
    return true;
  }
  else
    LOG_FAIL("replacement for throw call");
}
*/

Tunnel* Tunnel::GetSetTunnel()
{
  if (m_tunType == POSSTUNIN) {
    Tunnel *tun = (Tunnel*)GetNewInst();
    tun->m_tunType = SETTUNIN;
    return tun;
  }
  else if (m_tunType == POSSTUNOUT) {
    Tunnel *tun = (Tunnel*)GetNewInst();
    tun->m_tunType = SETTUNOUT;
    return tun;
  }
  else {
    cout << "GetSetTunnel on wrong tunnel type\n";
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

void Tunnel::Prop()
{
  if (!IsValidCost(m_cost)) {
#if DODLA
    DLANode::Prop();
#endif //DODLA

    if ((m_tunType == SETTUNIN || m_tunType == SETTUNOUT) && !m_pset) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    if ((m_tunType == POSSTUNIN && m_inputs.size() != 1)
	|| (m_tunType == POSSTUNOUT && m_children.size() != 1))
      {
	cout << "m_inputs.size() != 1\n";
	cout << "it's " << m_inputs.size() << endl;
	if (m_tunType == POSSTUNIN)
	  cout << "m_tunType == POSSTUNIN\n";
	else
	  cout << "m_tunType == POSSTUNOUT\n";
	LOG_FAIL("replacement for throw call");
	throw;
      }
    if (m_tunType == LASTTUNNEL) {
      cout << "bad tunnel type!\n";
      LOG_FAIL("replacement for throw call");
      throw;
    }
    if (m_tunType == SETTUNIN || m_tunType == POSSTUNIN)
      if (!m_inputs.size()) {
	cout << "!m_inputs.size() on this = " << this << "\n";
	if (m_tunType == POSSTUNIN)
	  cout << "m_tunType == POSSTUNIN\n";
	else 
	  cout << "m_tunType == SETTUNIN is on " << m_pset << endl;
	LOG_FAIL("replacement for throw call");
	throw;
      }


    if (m_tunType == POSSTUNIN) {
      if (m_inputs.size() != 1) {
	cout << "m_inputs.size() != 1\n";
	LOG_FAIL("replacement for throw call");
	throw;
      }
    }    
    else if (m_tunType == POSSTUNOUT) {
      if (!IsLoopTunnel() && m_inputs.size() != 1) {
	LOG_FAIL("replacement for throw call");
	throw;
      }
    }
    if (m_tunType == SETTUNIN || m_tunType == POSSTUNIN)
      if (!m_inputs.size()) {
	cout << "!m_inputs.size() on " << this << "\n";
	LOG_FAIL("replacement for throw call");
	throw;
      }

    m_cost = ZERO;
    if (m_tunType == SETTUNOUT) {
      m_pset->Prop();
    }
    else if (m_tunType == POSSTUNOUT || m_tunType == SETTUNIN) {
      for (ConnNum i = 0; i < m_inputs.size(); ++i) {
	Input(i)->Prop();
      }
    }
  }
}

const DataTypeInfo& Tunnel::DataType(ConnNum num) const
{
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (m_tunType == SETTUNOUT && !m_pset->IsReal()) {
    return GetRealTunnel()->DataType(num);
  }
  else {
    return InputDataType(0);
  }
}

#if TWOD
const SizeList* Tunnel::GetM(ConnNum num) const
{
  if (m_tunType == SETTUNOUT) {
    if (m_pset->IsReal())
      return ((DLANode*)Input(0))->GetM(num);
    else {
      return GetRealTunnel()->GetM(num);
    }
  }
  else if (m_tunType != POSSTUNIN) {
    return GetInputM(num);
  }
  else {
    return ((DLANode*)Input(0))->GetM(num);
  }
}

const SizeList* Tunnel::GetN(ConnNum num) const
{
  if (m_tunType == SETTUNOUT) {
    if (m_pset->IsReal())
      return ((DLANode*)Input(0))->GetN(num);
    else {
      return GetRealTunnel()->GetN(num);
    }
  }
  else if (m_tunType != POSSTUNIN) {
    return GetInputN(num);
  }
  else {
    return ((DLANode*)Input(0))->GetN(num);
  }
}

#if DODM
const SizeList* Tunnel::LocalM(ConnNum num) const
{
  if (m_tunType == SETTUNOUT && !m_pset->IsReal()) {
    return GetRealTunnel()->LocalM(num);
  }
  else {
    const NodeConn *conn = m_inputs[0];
    const DLANode *input = (DLANode*)(conn->m_n);
    
    return input->LocalM(conn->m_num);
  }
}

const SizeList* Tunnel::LocalN(ConnNum num) const
{
  if (m_tunType == SETTUNOUT && !m_pset->IsReal()) {
    return GetRealTunnel()->LocalN(num);
  }
  else {
    const NodeConn *conn = m_inputs[0];
    const DLANode *input = (DLANode*)(conn->m_n);
    
    return input->LocalN(conn->m_num);
  }
}
#endif

#elif DOTENSORS
const Dim Tunnel::NumDims(ConnNum num) const
{

  if (m_tunType == SETTUNOUT) {
    if (m_pset->IsReal())
      return ((DLANode*)Input(0))->NumDims(num);
    else {
      return GetRealTunnel()->NumDims(num);
    }
  }
  else if (m_tunType != POSSTUNIN) {
    return InputNumDims(num);
  }
  else {
    return ((DLANode*)Input(0))->InputNumDims(num);
  }
}

const SizeList* Tunnel::Len(ConnNum num,Dim dim) const
{

  if (m_tunType == SETTUNOUT) {
    if (m_pset->IsReal())
      return ((DLANode*)Input(0))->Len(num,dim);
    else {
      return GetRealTunnel()->Len(num,dim);
    }
  }
  else if (m_tunType != POSSTUNIN) {
    return InputLen(num,dim);
  }
  else {
    return ((DLANode*)Input(0))->Len(num,dim);
  }
}

const SizeList* Tunnel::LocalLen(ConnNum num,Dim dim) const
{
  if (m_tunType == SETTUNOUT && !m_pset->IsReal()) {
    return GetRealTunnel()->LocalLen(num,dim);
  }
  else {
    const NodeConn *conn = m_inputs[0];
    const DLANode *input = (DLANode*)(conn->m_n);
    
    return input->LocalLen(conn->m_num,dim);
  }
}

#endif

Name Tunnel::GetName(ConnNum num) const 
{
  if (m_tunType == SETTUNOUT) {
    if (m_pset->IsReal())
      return GetInputName(0);
    else {
      return GetRealTunnel()->GetName(num);
    }
  }
  else if (m_tunType != POSSTUNIN) {
    return GetInputName(num);
  }
  else {
    if (num > 0) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    return GetInputName(0);
  }
}

unsigned int Tunnel::NumOutputs() const
{
  switch(m_tunType) 
    {
    case (POSSTUNOUT) :
    case (SETTUNIN) :
      return m_inputs.size();
    case (POSSTUNIN) :
	return Input(0)->NumOutputs();
    case (SETTUNOUT) :
      if (m_pset->IsReal())
	return Input(0)->NumOutputs();
      else {
	return GetRealTunnel()->NumOutputs();
      }
    default:
      cout << "bad tunnel type\n";
      LOG_FAIL("replacement for throw call");
      throw;
      return -1;
  } 
}

void Tunnel::Duplicate(const Node *node, bool shallow, bool possMerging)
{
  TunnelSuper::Duplicate(node, shallow, possMerging);
  m_tunType = ((Tunnel*)node)->m_tunType;
}

string TunTypeToStr(TunType type)
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



void Tunnel::FlattenCore(ofstream &out) const
{
  TunnelSuper::FlattenCore(out);
  WRITE(m_tunType);
  WRITE(m_pset);
}



void Tunnel::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  TunnelSuper::UnflattenCore(in, info);
  READ(m_tunType);
  READ(m_pset);
  Swap(&m_pset,info.psetMap);
}

bool Tunnel::Overwrites(const Node *input, ConnNum num) const
{
  if (m_tunType == SETTUNIN) {
    if (m_pset->IsReal()) {
      if (!m_children.size()) {
	LOG_FAIL("replacement for throw call");
	throw;
      }
      const NodeConn *conn = m_children[0];
      return conn->m_n->Overwrites(this, 0);
    }
    else {
      const RealPSet *real = m_pset->GetReal();
      unsigned int num = FindInTunVec(m_pset->m_inTuns, this);
      const Node *realTun = real->m_inTuns[num];
      return realTun->Overwrites(NULL,0);
    }
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

Tunnel* Tunnel::GetRealTunnel() 
{
  if (!m_pset)
    return this;
  if (m_pset->IsReal())
    return this;
  RealPSet *real = m_pset->GetReal();
  if (m_tunType == SETTUNIN) {
    unsigned int num = FindInTunVec(m_pset->m_inTuns, this);
    return (Tunnel*)(real->m_inTuns[num]);
  }
  else if (m_tunType == SETTUNOUT) {
    unsigned int num = FindInTunVec(m_pset->m_outTuns, this);
    return (Tunnel*)(real->m_outTuns[num]);
  }
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
}


const Tunnel* Tunnel::GetRealTunnel() const
{
  if (!m_pset)
    return this;
  if (m_pset->IsReal())
    return this;
  const RealPSet *real = m_pset->GetReal();
  if (m_tunType == SETTUNIN) {
    unsigned int num = FindInTunVec(m_pset->m_inTuns, this);
    return (Tunnel*)(real->m_inTuns[num]);
  }
  else if (m_tunType == SETTUNOUT) {
    unsigned int num = FindInTunVec(m_pset->m_outTuns, this);
    return (Tunnel*)(real->m_outTuns[num]);
  }
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
}
