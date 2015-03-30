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



#include "combineSingleIter.h"
#include "splitSingleIter.h"
#include "elemRedist.h"
#include <cmath>

#if TWOD
CombineSingleIter::CombineSingleIter() 
  :CombineBase()
{
}

CombineSingleIter::CombineSingleIter(PartDir dir, TunType type) 
  : CombineBase(dir,type)
{
}
#else

CombineSingleIter::CombineSingleIter() 
  : CombineBase()
{
}

CombineSingleIter::CombineSingleIter(Dim partDim, TunType type) 
  : CombineBase(partDim, type)
{
}
#endif


void CombineSingleIter::Prop()
{
  if (!IsValidCost(m_cost)) {
    LoopTunnel::Prop();

#if TWOD
    if ((m_dir == PARTDOWN || m_dir == PARTUPWARD)
	&& (GetUpStat(TL) != GetUpStat(TR) || GetUpStat(BL) != GetUpStat(BR))) 
      {
	cout << "bad statuses\n";
	LOG_FAIL("replacement for throw call");
      }
    else if ((m_dir == PARTRIGHT || m_dir == PARTLEFT)
	     && (GetUpStat(TL) != GetUpStat(BL) || GetUpStat(TR) != GetUpStat(BR)))
      {
	cout << "bad statuses\n";
	LOG_FAIL("replacement for throw call");
      }
#endif
  
    if (m_tunType == SETTUNOUT) {
      if (!m_pset->IsLoop())
	LOG_FAIL("replacement for throw call");

      for(ConnNum i = 0; i < m_inputs.size(); ++i) {
	if (Input(i)->GetNodeClass() != CombineSingleIter::GetClass()) {
	  LOG_FAIL("replacement for throw call");
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
	  LOG_FAIL("replacement for throw call");
	}
        NodeConn *conn = m_inputs[m_inputs.size()-1];
#if TWOD
        if (conn->m_num != GetNumElems(m_dir)) {
#else
	  if (conn->m_num != 3) {
#endif
	    LOG_FAIL("replacement for throw call");
	  }
	  else if (conn->m_n->GetNodeClass() != SplitSingleIter::GetClass())
	    LOG_FAIL("replacement for throw call");
	
	if (Input(m_inputs.size()-1)->GetNodeClass() != SplitSingleIter::GetClass()) {
	  cout << "Last input isn't the right size\n";
	}
	for (unsigned int i = 0; i < m_children.size(); ++i) {
	  if (Child(i)->GetNodeClass() != CombineSingleIter::GetClass()) {
	    cout << Child(i)->GetType() << endl;
	    LOG_FAIL("replacement for throw call");
	  }
	}
      }
    else
	LOG_FAIL("replacement for throw call");

#if DODM
	const DataTypeInfo &type = GetRealTunnel()->InputDataType(0);
	for(ConnNum i = 0; i < m_inputs.size(); ++i) {
	  if (type != InputDataType(i)) {
	    cout << "Bad input types\n";
	    //      cout << InputDataType(i).Str() << endl;
	    //      cout << type.Str() << endl;
	    //From when we just used DistType w/out perm
	    //	    cout << DistTypeToStr(type) << endl;
	    //	    cout << DistTypeToStr(InputDataType(i).m_dist) << endl;
	    LOG_FAIL("replacement for throw call");
	  }
	}
#endif


    m_cost = ZERO;
  }
}
    
const DataTypeInfo& CombineSingleIter::DataType(ConnNum num) const
{
  if (m_tunType == POSSTUNOUT)
    return InputDataType(m_inputs.size()-1);
  else if (m_tunType == SETTUNOUT)
    return GetRealTunnel()->InputDataType(0);
  else
    LOG_FAIL("replacement for throw call");
}

#if TWOD
const Sizes* CombineSingleIter::GetM(ConnNum num) const
{
  if (num > 0)
    LOG_FAIL("replacement for throw call");
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(GetRealTunnel()->Input(0)->Input(GetNumElems(m_dir))->Input(0)))->GetInputM(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return GetInputM(m_inputs.size()-1);
  }
  else {
    return GetInputM(0);
  }
}

const Sizes* CombineSingleIter::GetN(ConnNum num) const
{
  if (num > 0)
    LOG_FAIL("replacement for throw call");
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(GetRealTunnel()->Input(0)->Input(GetNumElems(m_dir))->Input(0)))->GetInputN(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return GetInputN(m_inputs.size()-1);
  }
  else {
    return GetInputN(0);
  }
}

#if DODM
const Sizes* CombineSingleIter::LocalM(ConnNum num) const
{
  if (num > 0)
    LOG_FAIL("replacement for throw call");
  if (m_tunType == SETTUNOUT) {
    DLANode *possTunOut = (DLANode*)(GetRealTunnel()->Input(0));
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

const Sizes* CombineSingleIter::LocalN(ConnNum num) const
{
  if (num > 0)
    LOG_FAIL("replacement for throw call");
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(GetRealTunnel()->Input(0)->Input(GetNumElems(m_dir))->Input(0)))->InputLocalN(0);
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
const Dim CombineSingleIter::NumDims(ConnNum num) const
{
  if (num > 0)
    LOG_FAIL("replacement for throw call");
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(GetRealTunnel()->Input(0)->Input(3)->Input(0)))->NumDims(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputNumDims(m_inputs.size()-1);
  }
  else {
    return InputNumDims(0);
  }
}

const Sizes* CombineSingleIter::Len(ConnNum num, Dim dim) const
{
  if (num > 0)
    LOG_FAIL("replacement for throw call");
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(GetRealTunnel()->Input(0)->Input(3)->Input(0)))->Len(0,dim);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputLen(m_inputs.size()-1,dim);
  }
  else {
    return InputLen(0,dim);
  }
}

const Sizes* CombineSingleIter::LocalLen(ConnNum num, Dim dim) const
{
  if (num > 0)
    LOG_FAIL("replacement for throw call");
  if (m_tunType == SETTUNOUT) {
    DLANode *possTunOut = (DLANode*)(GetRealTunnel()->Input(0));
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

Name CombineSingleIter::GetName(ConnNum num) const
{
  if (num > 0)
    LOG_FAIL("replacement for throw call");
  if (m_tunType == POSSTUNOUT)
    return ((SplitSingleIter*)Input(m_inputs.size()-1))->GetOrigName();
  else
    return GetRealTunnel()->Input(0)->GetName(0);
}

Tunnel* CombineSingleIter::GetSetTunnel()
{
  CombineSingleIter *tun;
  /*
  if (m_tunType == POSSTUNIN)
#if TWOD
    tun = new CombineSingleIter(m_dir, SETTUNIN);
#else
    tun = new CombineSingleIter(m_partDim, SETTUNIN);
#endif
else*/
  if (m_tunType == POSSTUNOUT || m_tunType == SETTUNOUT)
#if TWOD
    tun = new CombineSingleIter(m_dir, SETTUNOUT);
#else
    tun = new CombineSingleIter(m_partDim, SETTUNOUT);
#endif
  else
    LOG_FAIL("replacement for throw call");
  tun->CopyTunnelInfo(this);
  return tun;
}

void CombineSingleIter::PrintCode(IndStream &out)
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
#elif DOTENSORS
  Name nameT = GetInputName(3);
  Name nameB = nameT;
  nameT.m_name += "_part" + std::to_string(m_partDim) + "T";
  nameB.m_name += "_part" + std::to_string(m_partDim) + "B";

      out.Indent();
      *out << "SlidePartitionDown\n"
	   << out.Tabs(0)
	   << "( " << nameT.str() << ",  " << GetInputNameStr(0) << ",\n"
	   << out.Tabs(0)
	   << "       " << GetInputNameStr(1) << ",\n"
	   << out.Tabs(0)
	   << "  /**/ /**/\n"
	   << out.Tabs(0)
	   << "  " << nameB.str() << ", " << GetInputNameStr(2) << ", " 
	   << m_partDim << " );\n";
#else
  LOG_FAIL("replacement for throw call");
#endif
}

NodeType CombineSingleIter::GetType() const
{
#if TWOD
  return "Combine " + PartDirToStr(m_dir) + "( " + Tunnel::GetType() + " )";
#else
  string str = "Combine ";
  str += m_partDim;
  return str + " ( " + Tunnel::GetType() + " )";
#endif
}


LoopTunnel* CombineSingleIter::GetMatchingInTun() const
{
  if (m_tunType == SETTUNOUT)
    if (m_pset->IsReal())
      return (LoopTunnel*)(((LoopTunnel*)Input(0))->GetMatchingInTun()->Input(0));
    else {
      LoopTunnel *possTunIn = ((LoopTunnel*)(GetRealTunnel()->Input(0)))->GetMatchingInTun();
      LoopTunnel *realSetTunIn = (LoopTunnel*)(possTunIn->Input(0));
      return (LoopTunnel*)(m_pset->m_inTuns[FindInTunVec(m_pset->GetReal()->m_inTuns, realSetTunIn)]);
    }
  else if (m_tunType != POSSTUNOUT)
    LOG_FAIL("replacement for throw call");
  
#if DOTENSORS
  const Node *in = Input(3);
#else
  const Node *in = Input(GetNumElems(m_dir));
#endif
  if (in->IsTunnel(POSSTUNIN) && ((Tunnel*)in)->IsSplit()) {
    return (LoopTunnel*)in;
  }
  else {
    cout << "Didn't find matching in tun\n";
    LOG_FAIL("replacement for throw call");
  }
}
