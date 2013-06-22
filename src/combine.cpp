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



#include "combine.h"
#include "split.h"
#include "distributions.h"
#include <cmath>

void Combine::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (m_tunType==POSSTUNOUT) {
      if (m_inputs.size() != GetNumElems(m_dir)+1) {
        cout << "wrong number of inputs " << this << endl;
        cout << m_inputs.size() << " inputs instead of " << GetNumElems(m_dir)+1 << endl;
        throw;
      }
      else {
        NodeConn *conn = m_inputs[m_inputs.size()-1];
        if (conn->m_num != GetNumElems(m_dir))
          throw;
        else if (conn->m_n->GetNodeClass() != Split::GetClass())
          throw;
      }
    }
    
    if (m_tunType==POSSTUNOUT || GetMyLoop()->GetType()==ELEMLOOP) {
      unsigned int size = m_inputs.size();
      for (unsigned int i = 0; i < size; ++i) {
        Node *input = Input(i);
        input->Prop();
      }
    }
    
    m_cost = ZERO;
    DistType type = InputDistType(0);
    for(unsigned int i = 0; i < m_inputs.size(); ++i) {
      if (type != InputDistType(i)) {
        cout << "Bad input types\n";
        cout << DistTypeToStr(type) << endl;
        cout << DistTypeToStr(InputDistType(i)) << endl;
        throw;
        m_poss->MarkInsane();
      }
    }
  }
}

const Sizes* Combine::GetM(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(Input(0)->Input(GetNumElems(m_dir))->Input(0)))->GetInputM(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return GetInputM(m_inputs.size()-1);
  }
  else {
    return GetInputM(0);
  }
}

const Sizes* Combine::GetN(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(Input(0)->Input(GetNumElems(m_dir))->Input(0)))->GetInputN(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return GetInputN(m_inputs.size()-1);
  }
  else {
    return GetInputN(0);
  }
}

const Sizes* Combine::LocalM(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    DLANode *possTunOut = (DLANode*)Input(0);
    DLANode *possTunIn = (DLANode*)(possTunOut->Input(GetNumElems(m_dir)));
    DLANode *setTunIn = (DLANode*)(possTunIn->Input(0));
    return setTunIn->InputLocalM(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputLocalM(m_inputs.size()-1);
  }
  else {
    return InputLocalM(0);
  }
}

const Sizes* Combine::LocalN(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(Input(0)->Input(GetNumElems(m_dir))->Input(0)))->InputLocalN(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputLocalN(m_inputs.size()-1);
  }
  else {
    return InputLocalN(0);
  }
}

Name Combine::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == POSSTUNOUT)
    return ((Split*)Input(m_inputs.size()-1))->GetOrigName();
  else
    return Input(0)->GetName(0);
}

PossTunnel* Combine::GetSetTunnel()
{
  Combine *tun;
  if (m_tunType == POSSTUNIN)
    tun = new Combine(m_dir, SETTUNIN);
  else if (m_tunType == POSSTUNOUT)
    tun = new Combine(m_dir, SETTUNOUT);
  else
    throw;
  tun->CopyTunnelInfo(this);
  return tun;
}

void Combine::PrintCode(IndStream &out)
{
  if (m_tunType != POSSTUNOUT) {
    //    cout << "returning from " << GetNameStr(0) << endl;
    return;
  }
  LoopType loopType = GetLoopType();
  if (loopType == ELEMLOOP) {
    switch(m_dir) {
    case (PARTDOWN):
      out.Indent();
      *out << "SlidePartitionDown\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(3) << "T,  " << GetInputNameStr(0) << ",\n"
	   << out.Tabs(0)
	   << "       " << GetInputNameStr(1) << ",\n"
	   << out.Tabs(0)
	   << "  /**/ /**/\n"
	   << out.Tabs(0)
	   << "  " << GetInputNameStr(3) << "B, " << GetInputNameStr(2) << " );\n";
      break;
    case (PARTUPWARD):
      out.Indent();
      *out << "SlidePartitionUp\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(3) << "T,  " << GetInputNameStr(0) << ",\n"
	   << out.Tabs(0)
	   << "  /**/ /**/\n"
	   << out.Tabs(0)
	   << "       " << GetInputNameStr(1) << ",\n"
	   << out.Tabs(0)
	   << "  " << GetInputNameStr(3) << "B, " << GetInputNameStr(2) << " );\n";
      break;
    case (PARTRIGHT):
      out.Indent();
      *out << "SlidePartitionRight\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(3) << "L,      " << GetInputNameStr(3) << "R,\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(0) << ", " << GetInputNameStr(1) << ", /**/ " << GetInputNameStr(2) << " );\n";
      break;
    case (PARTLEFT):
      out.Indent();
      *out << "SlidePartitionLeft\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(3) << "L,      " << GetInputNameStr(3) << "R,\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(0) << ", /**/ " << GetInputNameStr(1) << ", " << GetInputNameStr(2) << " );\n";
      break;
    case (PARTDIAG):
      out.Indent();
      *out << "SlidePartitionDownDiagonal\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(9) << "TL, /**/ " << GetInputNameStr(9) << "TR,  "
	   << GetInputNameStr(0) << ", " << GetInputNameStr(3) << ", /**/ " << GetInputNameStr(6) << ",\n"
	   << out.Tabs(0)
	   << "       /**/       "
	   << GetInputNameStr(1) << ", " << GetInputNameStr(4) << ", /**/ " << GetInputNameStr(7) << ",\n"
	   << out.Tabs(0)
	   << " /*************/ /******************/\n"
	   << out.Tabs(0)
	   << "  " << GetInputNameStr(9) << "BL, /**/ " << GetInputNameStr(9) << "BR,  "
	   << GetInputNameStr(2) << ", " << GetInputNameStr(5) << ", /**/ " << GetInputNameStr(8) << " );\n";
      break;
    case (PARTDIAGBACK):
      out.Indent();
      *out << "SlidePartitionUpDiagonal\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(9) << "TL, /**/ " << GetInputNameStr(9) << "TR,  "
	   << GetInputNameStr(0) << ", /**/ " << GetInputNameStr(3) << ", " << GetInputNameStr(6) << ",\n"
	   << out.Tabs(0)
	   << " /*************/ /******************/\n"
	   << out.Tabs(0)
	   << GetInputNameStr(1) << ", /**/ " << GetInputNameStr(4) << ", " << GetInputNameStr(7) << ",\n"
	   << "       /**/       "
	   << out.Tabs(0)
	   << "  " << GetInputNameStr(9) << "BL, /**/ " << GetInputNameStr(9) << "BR,  "
	   << GetInputNameStr(2) << ", /**/ " << GetInputNameStr(5) << ", " << GetInputNameStr(8) << " );\n";
      break;
    default:
      cout << "bad dir\n";
      break;
    }
  }
}

void Combine::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  LoopTunnel::Duplicate(orig, shallow, possMerging);
  const Combine *com = (Combine*)orig;
  m_dir = com->m_dir;
}

NodeType Combine::GetType() const
{
  return "Combine " + PartDirToStr(m_dir) + "( " + PossTunnel::GetType() + " )";
}

void Combine::SanityCheck()
{
  LoopTunnel::SanityCheck();
  
  if ((m_dir == PARTDOWN || m_dir == PARTUPWARD)
      && (GetUpStat(TL) != GetUpStat(TR) || GetUpStat(BL) != GetUpStat(BR))) 
  {
    cout << "bad statuses\n";
    throw;
  }
  else if ((m_dir == PARTRIGHT || m_dir == PARTLEFT)
           && (GetUpStat(TL) != GetUpStat(BL) || GetUpStat(TR) != GetUpStat(BR)))
  {
    cout << "bad statuses\n";
    throw;
  }
  
  if (m_tunType == SETTUNOUT) {
    if (!m_pset->IsLoop())
      throw;
    
    for(unsigned int i = 0; i < m_inputs.size(); ++i) {
      if (Input(i)->GetNodeClass() != Combine::GetClass()) {
        throw;
      }
    }
  }
  else if (m_tunType == POSSTUNOUT) {
    if (m_inputs.size() != GetNumElems(m_dir)+1) {
      int num = m_inputs.size();
      cout << "combine has wrong number of inputs\n";
      cout << num <<endl;
      throw;
    }
    if (Input(m_inputs.size()-1)->GetNodeClass() != Split::GetClass()) {
      cout << "Last input isn't the right size\n";
    }
    for (unsigned int i = 0; i < m_children.size(); ++i) {
      if (Child(i)->GetNodeClass() != Combine::GetClass()) {
        cout << Child(i)->GetType() << endl;
        throw;
      }
    }
  }
  else
    throw;
}

void Combine::FlattenCore(ofstream &out) const
{
  LoopTunnel::FlattenCore(out);
  WRITE(m_dir);
}

void Combine::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  LoopTunnel::UnflattenCore(in,info);
  READ(m_dir);
}
