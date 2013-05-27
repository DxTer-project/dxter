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



#include "loopTunnel.h"
#include "distributions.h"
#include <cmath>


LoopTunnel::LoopTunnel(PossTunType type)
: PossTunnel(type)
{
  m_statTL = BADUP;
  m_statTR = BADUP;
  m_statBL = BADUP;
  m_statBR = BADUP;
  m_msizes = NULL;
  m_nsizes = NULL;
  m_mlsizes = NULL;
  m_nlsizes = NULL;
}

LoopTunnel::~LoopTunnel()
{
  if (m_msizes) {
    delete m_msizes;
    m_msizes = NULL;
    delete m_mlsizes;
    m_mlsizes = NULL;
    delete m_nsizes;
    m_nsizes = NULL;
    delete m_nlsizes;
    m_nlsizes = NULL;
  }
}

void LoopTunnel::SetUpStats( UpStat statTL, UpStat statTR,
                            UpStat statBL, UpStat statBR )
{
  m_statTL = statTL;
  m_statTR = statTR;
  m_statBL = statBL;
  m_statBR = statBR;
}

void LoopTunnel::SetAllStats(UpStat stat)
{
  SetUpStats(stat,stat,stat,stat);
}

void LoopTunnel::CopyUpStats(const LoopTunnel *tun)
{
  m_statTL = tun->m_statTL;
  m_statTR = tun->m_statTR;
  m_statBL = tun->m_statBL;
  m_statBR = tun->m_statBR;
}

void LoopTunnel::SanityCheck()
{
  PossTunnel::SanityCheck();
  if (m_statTL == BADUP || m_statTR == BADUP
      || m_statBL == BADUP || m_statBR == BADUP)
  {
    cout << "update status not set for " << GetType() << " " << this << endl;
    cout << "it's a " << TunTypeToStr(m_tunType) << endl;
    throw;
  }
  if (m_tunType == POSSTUNOUT) {
    Node *in = m_inputs[m_inputs.size()-1]->m_n;
    if (!in->IsLoopTunnel())
      throw;
    if (((LoopTunnel*)in)->m_tunType != POSSTUNIN)
      throw;
  }
}

PossTunnel* LoopTunnel::GetSetTunnel()
{
  LoopTunnel *tun;
  if (m_tunType == POSSTUNIN)
    tun = new LoopTunnel(SETTUNIN);
  else if (m_tunType == POSSTUNOUT)
    tun = new LoopTunnel(SETTUNOUT);
  else
    throw;
  tun->CopyUpStats(this);
  return tun;
}

void LoopTunnel::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (m_tunType == SETTUNIN && m_inputs.size() != 1) {
      cout << "1 m_inputs.size() = " << m_inputs.size() << endl;
      for(unsigned int i = 0; i < m_inputs.size(); ++i) {
        cout << Input(i)->GetType() << " " << Input(i)->GetNameStr(0) << endl;
      }
      throw;
    }
    if (m_tunType == POSSTUNOUT && m_inputs.size() != 2) {
      if (Input(1)->GetNodeClass() != LoopTunnel::GetClass())
        throw;
      cout << "2 m_inputs.size() = " << m_inputs.size() << endl;
      throw;
      Input(1)->Prop();
    }
    else if (m_tunType == POSSTUNIN) {
      if (Input(0)->GetNodeClass() != LoopTunnel::GetClass()) {
        cout << Input(0)->GetNodeClass() << endl;
        throw;
      }
    }
    m_cost = ZERO;
    if (m_tunType==POSSTUNOUT || GetMyLoop()->GetType()==ELEMLOOP) {
      if (m_inputs.size()) 
        Input(0)->Prop();
    }
  }
}

const Sizes* LoopTunnel::GetM(unsigned int num) const
{
  switch(m_tunType) 
  {
      //    case (SETTUNIN):
    case (POSSTUNOUT):
      if (num > 0)
        throw;
      return ((DLANode*)(Input(1)->Input(0)))->GetInputM(0);
    case (SETTUNOUT):
      if (num > 0)
        throw;
      return GetInputM(0);
    case (POSSTUNIN):
      if (num == 0) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        return input->m_msizes;
      }
      else if (num == 1) {
        return GetInputM(0);
      }
      else
        throw;
    default:
      throw;
  }
}

const Sizes* LoopTunnel::GetN(unsigned int num) const
{
  switch(m_tunType) 
  {
    case (SETTUNIN):
    case (POSSTUNOUT):
      if (num > 0)
        throw;
      return ((DLANode*)(Input(1)->Input(0)))->GetInputN(0);
    case (SETTUNOUT):
      if (num > 0)
        throw;
      return GetInputN(0);
    case (POSSTUNIN):
      if (num == 0) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        return input->m_nsizes;
      }
      else if (num == 1) {
        return GetInputN(0);
      }
      else
        throw;
    default:
      throw;
  }
}

const Sizes* LoopTunnel::LocalM(unsigned int num) const
{
  switch(m_tunType) 
  {
    case (SETTUNIN):
    case (POSSTUNOUT):
      if (num > 0)
        throw;
      return ((DLANode*)(Input(1)->Input(0)))->InputLocalM(0);
    case (SETTUNOUT):
      if (num > 0)
        throw;
      return InputLocalM(0);
    case (POSSTUNIN):
      if (num == 0) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        return input->m_mlsizes;
      }
      else if (num == 1) {
        return InputLocalM(0);
      }
      else
        throw;
    default:
      throw;
  }
}

const Sizes* LoopTunnel::LocalN(unsigned int num) const
{
  switch(m_tunType) 
  {
    case (SETTUNIN):
    case (POSSTUNOUT):
      if (num > 0)
        throw;
      return ((DLANode*)(Input(1)->Input(0)))->InputLocalN(0);
    case (SETTUNOUT):
      if (num > 0)
        throw;
      return InputLocalN(0);
    case (POSSTUNIN):
      if (num == 0) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        return input->m_nlsizes;
      }
      else if (num == 1) {
        return InputLocalN(0);
      }
      else
        throw;
    default:
      throw;
  }
}

unsigned int LoopTunnel::NumOutputs() const
{
  switch(m_tunType) 
  {
    case (POSSTUNOUT) :
    case (SETTUNIN) :
    case (SETTUNOUT) :
      return 1;
    case (POSSTUNIN) :
      return 2;
    default:
      cout << "bad tunnel type\n";
      throw;
      return -1;
  } 
}

Name LoopTunnel::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  Name name;
  if (GetLoopType() == BLISLOOP) {
    if (m_tunType == SETTUNOUT)
      return ((LoopTunnel*)Input(0))->GetOrigName();
    else if (m_tunType == POSSTUNOUT)
      return ((LoopTunnel*)Input(1))->GetOrigName();
    else if (m_tunType == SETTUNIN) {
      name.m_name = GetInputName(0).str();
//      if (name.m_name[name.m_name.length()-1] != '_')
//        name.m_name = name.m_name + "_";
      return name;
    }
  }
  return GetInputName(0);
}

Name LoopTunnel::GetOrigName() const
{
  if (m_tunType == SETTUNOUT || m_tunType == POSSTUNIN)
    return ((LoopTunnel*)Input(0))->GetOrigName();
  else if (m_tunType == POSSTUNOUT)
    return ((LoopTunnel*)Input(m_inputs.size()-1))->GetOrigName();
  else if (m_tunType == SETTUNIN)
    return ((LoopTunnel*)Input(0))->GetName(InputConnNum(0));
  else
    throw;
}

Loop* LoopTunnel::GetMyLoop() const
{
  PSet *set = NULL;
  if (m_tunType == SETTUNIN)
    set = m_pset;
  else if (m_tunType == SETTUNOUT)
    set = m_pset;
  else if (m_tunType == POSSTUNIN)
    set = ((PossTunnel*)Input(0))->m_pset;
  else if (m_tunType == POSSTUNOUT)
    set = ((PossTunnel*)Child(0))->m_pset;
  else 
    throw;
  if (!set->IsLoop()) {
    cout << "Loop Tunnel doesn't have Loop as PSet\n";
    throw;
  }
  return (Loop*)set;
}


void LoopTunnel::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  PossTunnel::Duplicate(orig, shallow, possMerging);
  CopyUpStats((LoopTunnel*)orig);
}

NodeType LoopTunnel::GetType() const
{
  return "LoopTunnel (" + PossTunnel::GetType() + ")";
}

UpStat LoopTunnel::GetUpStat(Quad quad) const
{
  if (m_tunType == SETTUNIN || m_tunType == SETTUNOUT) {
    switch(quad)
    {
      case(TL):
        return m_statTL;
      case(TR):
        return m_statTR;
      case(BL):
        return m_statBL;
      case(BR):
        return m_statBR;
      default:
        throw;
    }
  }
  else if (m_tunType == POSSTUNIN) {
    return ((LoopTunnel*)Input(0))->GetUpStat(quad);
  }
  else if (m_tunType == POSSTUNOUT) {
    return ((LoopTunnel*)Child(0))->GetUpStat(quad);
  }
  else
    throw;
}

bool LoopTunnel::AllFullyUpdated() const
{
  return (m_statTL == FULLUP &&
          m_statTR == FULLUP &&
          m_statBL == FULLUP &&
          m_statBR == FULLUP);
}

bool LoopTunnel::QuadInUse(Quad quad, bool atEnd) const
{
  if (m_tunType == SETTUNIN) {
    Node *child = Child(0);
    if (!child->IsLoopTunnel())
      throw;
    return ((LoopTunnel*)child)->QuadInUse(quad, atEnd);
  }
  else if (m_tunType == POSSTUNIN) {
    NodeConnVecConstIter iter = m_children.begin();
    for( ; iter != m_children.end(); ++iter) {
      if (!(*iter)->m_n->IsPossTunnel(POSSTUNOUT))
        return true;
    }
    return false;
  }
  else
    throw;
}

LoopTunnel* LoopTunnel::GetMatchingOutTun() const
{
  if (m_tunType == SETTUNIN)
    return (LoopTunnel*)(((LoopTunnel*)Child(0))->GetMatchingOutTun()->Child(0));
  else if (m_tunType != POSSTUNIN) 
    throw;
  
  NodeConnVecConstIter iter = m_children.begin();
  for(; iter != m_children.end(); ++iter) {
    const NodeConn *con = *iter;
    const Node *child = con->m_n;
    if (child->IsPossTunnel(POSSTUNOUT) &&
        child->IsLoopTunnel()) {
      return (LoopTunnel*)child;
    }
  }
  cout << "Didn't find matching out tun\n";
  throw;
}

LoopTunnel* LoopTunnel::GetMatchingInTun() const
{
  if (m_tunType == SETTUNOUT)
    return (LoopTunnel*)(((LoopTunnel*)Input(0))->GetMatchingInTun()->Input(0));
  else if (m_tunType != POSSTUNOUT)
    throw;
  
  const Node *in = Input(1);
  if (in->IsPossTunnel(POSSTUNIN) && in->IsLoopTunnel()) {
    return (LoopTunnel*)in;
  }
  else {
    cout << "Didn't find matching in tun\n";
    throw;
  }
}

bool LoopTunnel::IsConst() const
{
  return AllFullyUpdated();
}

LoopType LoopTunnel::GetLoopType() const
{
  if (m_tunType == POSSTUNIN) {
    const Node *in = Input(0);
      if (!in->IsLoopTunnel())
        throw;
    return ((LoopTunnel*)in)->GetLoopType();
  }
  else if (m_tunType == POSSTUNOUT) {
    const Node *child = Child(0);
    if (!child->IsLoopTunnel())
      throw;
    return ((LoopTunnel*)child)->GetLoopType();
  }
  if (!m_pset) {
    throw;
  }
  return ((Loop*)m_pset)->GetType();
}

void LoopTunnel::FlattenCore(ofstream &out) const
{
  PossTunnel::FlattenCore(out);
  WRITE(m_statTL);
  WRITE(m_statTR);
  WRITE(m_statBL);
  WRITE(m_statBR);
}

void LoopTunnel::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  PossTunnel::UnflattenCore(in,info);
  READ(m_statTL);
  READ(m_statTR);
  READ(m_statBL);
  READ(m_statBR);
}

void LoopTunnel::StartFillingSizes()
{
  if (m_msizes)
    throw;
  if (m_tunType != SETTUNIN)
    return;
  m_msizes = new Sizes;
  m_mlsizes = new Sizes;
  m_nsizes = new Sizes;
  m_nlsizes = new Sizes;
}

void LoopTunnel::AppendSizes(unsigned int execNum, unsigned int numIters)
{
  if (m_tunType != SETTUNIN)
    return;
  if (!m_msizes)
    throw;
  const DLANode *input = (DLANode*)Input(0);
  unsigned int num = InputConnNum(0);
  const Sizes *ms = input->GetM(num);
  if (!ms) {
    cout << "!ms\n";
    throw;
  }
  const Sizes *ns = input->GetN(num);
  const Sizes *lms = input->LocalM(num);
  const Sizes *lns = input->LocalN(num);
  unsigned int length = ms->NumSizes();
  if (length != ns->NumSizes() 
      || length != lms->NumSizes() 
      || length != lns->NumSizes() 
      || length <= execNum) 
  {
    cout << ms->NumSizes() << endl;
    cout << ns->NumSizes() << endl;
    cout << lms->NumSizes() << endl;
    cout << lns->NumSizes() << endl;
    throw;
  }
  const Size m = (*ms)[execNum];
  const Size n = (*ns)[execNum];
  const Size lm = (*lms)[execNum];
  const Size ln = (*lns)[execNum];
  m_msizes->AddRepeatedSizes(m, numIters);
  m_nsizes->AddRepeatedSizes(n, numIters);
  m_mlsizes->AddRepeatedSizes(lm, numIters);
  m_nlsizes->AddRepeatedSizes(ln, numIters);
}

void LoopTunnel::ClearSizeCache()
{
  if (!m_msizes) 
    return;
  delete m_msizes;
  m_msizes = NULL;
  delete m_mlsizes;
  m_mlsizes = NULL;
  delete m_nsizes;
  m_nsizes = NULL;
  delete m_nlsizes;
  m_nlsizes = NULL;
}

bool LoopTunnel::Overwrites(const Node *input, unsigned int num) const
{
  if (m_tunType == POSSTUNIN || m_tunType == SETTUNIN)
    return !IsConst();
  else
    return false;
}

bool LoopTunnel::KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const
{
  if (Input(m_inputs.size()-1) == input && InputConnNum(m_inputs.size()-1) == numIn) {
    numOut = numIn;
    return true;
  }
  else {
    return false;
  }
}
