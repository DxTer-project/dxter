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



#include "combineUnrolled.h"
#include "splitUnrolled.h"
#include "elemRedist.h"
#include <cmath>

#if TWOD
CombineUnrolled::CombineUnrolled() 
  : CombineBase (),
    m_unrollFactor(0)
{
}

CombineUnrolled::CombineUnrolled(PartDir dir, unsigned int unrollFactor, PossTunType type) 
  : CombineBase(dir,type), 
    m_unrollFactor(unrollFactor)
{
}

#else
CombineUnrolled::CombineUnrolled() 
  : CombineBase(),
   m_unrollFactor(0)
{
}


CombineUnrolled::CombineUnrolled(Dim partDim, unsigned int unrollFactor, PossTunType type) 
  : CombineBase(partDim, type), 
    m_unrollFactor(unrollFactor)
{
}
#endif


void CombineUnrolled::Prop()
{
  if (!IsValidCost(m_cost)) {
    CombineBase::Prop();
    
    if (m_tunType == POSSTUNOUT) {
      if (m_inputs.size() != m_unrollFactor+1)
	throw;
      
      if (Input(m_unrollFactor)->GetNodeClass() != SplitUnrolled::GetClass())
	throw;
    }

    if (!m_unrollFactor)
      throw;

    m_cost = ZERO;
  }
}

const DataTypeInfo& CombineUnrolled::DataType(unsigned int num) const
{
  return InputDataType(0);
}

#if TWOD
const Sizes* CombineUnrolled::GetM(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(Input(0)->Input(m_unrollFactor)->Input(0)))->GetInputM(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return GetInputM(m_unrollFactor);
  }
  else {
    return GetInputM(0);
  }
}

const Sizes* CombineUnrolled::GetN(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(Input(0)->Input(m_unrollFactor)->Input(0)))->GetInputN(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return GetInputN(m_unrollFactor);
  }
  else {
    return GetInputN(0);
  }
}

#if DODM
const Sizes* CombineUnrolled::LocalM(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    DLANode *possTunOut = (DLANode*)Input(0);
    DLANode *possTunIn = (DLANode*)(possTunOut->Input(m_unrollFactor));
    DLANode *setTunIn = (DLANode*)(possTunIn->Input(0));
    return setTunIn->InputLocalM(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputLocalM(m_unrollFactor);
  }
  else {
    return InputLocalM(0);
  }
}

const Sizes* CombineUnrolled::LocalN(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(Input(0)->Input(m_unrollFactor)->Input(0)))->InputLocalN(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputLocalN(m_unrollFactor);
  }
  else {
    return InputLocalN(0);
  }
}
#endif //DODM



#else
const Dim CombineUnrolled::NumDims(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(Input(0)->Input(m_unrollFactor)->Input(0)))->NumDims(0);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputNumDims(m_inputs.size()-1);
  }
  else {
    return InputNumDims(0);
  }
}

const Sizes* CombineUnrolled::Len(unsigned int num, Dim dim) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    return ((DLANode*)(Input(0)->Input(m_unrollFactor)->Input(0)))->Len(0,dim);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputLen(m_unrollFactor,dim);
  }
  else {
    return InputLen(0,dim);
  }
}

const Sizes* CombineUnrolled::LocalLen(unsigned int num, Dim dim) const
{
  if (num > 0)
    throw;
  if (m_tunType == SETTUNOUT) {
    DLANode *possTunOut = (DLANode*)Input(0);
    DLANode *possTunIn = (DLANode*)(possTunOut->Input(m_unrollFactor));
    DLANode *setTunIn = (DLANode*)(possTunIn->Input(0));
    return setTunIn->InputLocalLen(0,dim);
  }
  else if (m_tunType == POSSTUNOUT) {
    return InputLocalLen(m_unrollFactor,dim);
  }
  else {
    return InputLocalLen(0,dim);
  }
}


#endif

Name CombineUnrolled::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_tunType == POSSTUNOUT)
    return ((SplitUnrolled*)Input(m_unrollFactor))->GetOrigName();
  else
    return Input(0)->GetName(0);
}

PossTunnel* CombineUnrolled::GetSetTunnel()
{
  CombineUnrolled *tun;
  if (m_tunType == POSSTUNIN)
#if TWOD
    tun = new CombineUnrolled(m_dir, m_unrollFactor, SETTUNIN);
#else
    tun = new CombineUnrolled(m_partDim, m_unrollFactor, SETTUNIN);
#endif
  else if (m_tunType == POSSTUNOUT)
#if TWOD
    tun = new CombineUnrolled(m_dir, m_unrollFactor, SETTUNOUT);
#else
    tun = new CombineUnrolled(m_partDim, m_unrollFactor, SETTUNOUT);
#endif
  else
    throw;
  tun->CopyTunnelInfo(this);
  return tun;
}

void CombineUnrolled::PrintCode(IndStream &out)
{
  if (m_tunType != POSSTUNOUT) {
    //    cout << "returning from " << GetNameStr(0) << endl;
    return;
  }
#if TWOD
  LoopType loopType = GetLoopType();
  if (loopType == ELEMLOOP) {
    throw;
  }
#else
  throw;
#endif
}

void CombineUnrolled::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  CombineBase::Duplicate(orig, shallow, possMerging);
  const CombineUnrolled *com = (CombineUnrolled*)orig;
  m_unrollFactor = com->m_unrollFactor;
}

NodeType CombineUnrolled::GetType() const
{
#if TWOD
  return "CombineUnroled " + std::to_string(m_unrollFactor) 
    + PartDirToStr(m_dir) + "( " + PossTunnel::GetType() + " )";
#else
  string str = "CombineUnroled " + std::to_string(m_unrollFactor);
  str += m_partDim;
  return str + " ( " + PossTunnel::GetType() + " )";
#endif
}

void CombineUnrolled::FlattenCore(ofstream &out) const
{
  CombineBase::FlattenCore(out);
  WRITE(m_unrollFactor);
}

void CombineUnrolled::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  CombineBase::UnflattenCore(in,info);
  READ(m_unrollFactor);
}
