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



#include "combine.h"
#include "splitSingleIter.h"
#include "elemRedist.h"
#include <cmath>

void Combine::Prop()
{
  if (!IsValidCost(m_cost)) {
    LoopTunnel::Prop();

#if TWOD
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
#else
    throw;
#endif
  
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
#if TWOD
      if (m_inputs.size() != GetNumElems(m_dir)+1) {
#else
	if (m_inputs.size() != 4) {
#endif
	  int num = m_inputs.size();
	  cout << "combine has wrong number of inputs\n";
	  cout << num <<endl;
	  throw;
	}
        NodeConn *conn = m_inputs[m_inputs.size()-1];
#if TWOD
        if (conn->m_num != GetNumElems(m_dir)) {
#else
	  if (conn->m_num != 3) {
#endif
	    throw;
	  }
	  else if (conn->m_n->GetNodeClass() != SplitSingleIter::GetClass())
	    throw;
	
	if (Input(m_inputs.size()-1)->GetNodeClass() != SplitSingleIter::GetClass()) {
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

#if DODM
	DistType type = this->InputDataType(0).m_dist;
    for(unsigned int i = 0; i < m_inputs.size(); ++i) {
      if (DistTypeNotEqual(type, InputDataType(i).m_dist)) {
        cout << "Bad input types\n";
        cout << DistTypeToStr(type) << endl;
        cout << DistTypeToStr(InputDataType(i).m_dist) << endl;
        throw;
      }
    }
#endif


    m_cost = ZERO;
  }
}

    const DataTypeInfo& Combine::DataType(unsigned int num) const
    {
      return InputDataType(0);
    }

#if TWOD
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

#if DODM
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
#endif //DODM



#else
const Dim Combine::NumDims(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(Input(0)->Input(3)->Input(0)))->NumDims(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputNumDims(m_inputs.size()-1);
  }
  else {
    return InputNumDims(0);
  }
}

const Sizes* Combine::Len(unsigned int num, Dim dim) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(Input(0)->Input(3)->Input(0)))->Len(0,dim);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputLen(m_inputs.size()-1,dim);
  }
  else {
    return InputLen(0,dim);
  }
}

const Sizes* Combine::LocalLen(unsigned int num, Dim dim) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    DLANode *possTunOut = (DLANode*)Input(0);
#if TWOD
    DLANode *possTunIn = (DLANode*)(possTunOut->Input(GetNumElems(m_dir)));
#else
    DLANode *possTunIn = (DLANode*)(possTunOut->Input(3));
#endif
    DLANode *setTunIn = (DLANode*)(possTunIn->Input(0));
    return setTunIn->InputLocalLen(0,dim);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputLocalLen(m_inputs.size()-1,dim);
  }
  else {
    return InputLocalLen(0,dim);
  }
}


#endif

Name Combine::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == POSSTUNOUT)
    return ((SplitSingleIter*)Input(m_inputs.size()-1))->GetOrigName();
  else
    return Input(0)->GetName(0);
}

PossTunnel* Combine::GetSetTunnel()
{
  Combine *tun;
  if (m_tunType == POSSTUNIN)
#if TWOD
    tun = new Combine(m_dir, SETTUNIN);
#else
    tun = new Combine(m_partDim, SETTUNIN);
#endif
  else if (m_tunType == POSSTUNOUT)
#if TWOD
    tun = new Combine(m_dir, SETTUNOUT);
#else
    tun = new Combine(m_partDim, SETTUNOUT);
#endif
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
#if TWOD
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
#else
  throw;
#endif
}

void Combine::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  LoopTunnel::Duplicate(orig, shallow, possMerging);
  const Combine *com = (Combine*)orig;
#if TWOD
  m_dir = com->m_dir;
#else
  m_partDim = com->m_partDim;
#endif
}

NodeType Combine::GetType() const
{
#if TWOD
  return "Combine " + PartDirToStr(m_dir) + "( " + PossTunnel::GetType() + " )";
#else
  string str = "Combine ";
  str += m_partDim;
  return str + " ( " + PossTunnel::GetType() + " )";
#endif
}

void Combine::FlattenCore(ofstream &out) const
{
  LoopTunnel::FlattenCore(out);
#if TWOD
  WRITE(m_dir);
#else
  WRITE(m_partDim);
#endif
}

void Combine::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  LoopTunnel::UnflattenCore(in,info);
#if TWOD
  READ(m_dir);
#else
  READ(m_partDim);
#endif
}
