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

#include "splitSingleIter.h"
#include "combineSingleIter.h"
#include "elemRedist.h"
#include <cmath>
#include "LLDLA.h"
#include "logging.h"

#include "helperNodes.h"

#include "partition.h"

#if TWOD
SplitSingleIter::SplitSingleIter() 
  : SplitBase()
{
  m_addDir = false;
}

SplitSingleIter::SplitSingleIter(PartDir dir, TunType type, bool isControl) 
  : SplitBase(dir, type, isControl)
{
  m_addDir = false;
}
#else

SplitSingleIter::SplitSingleIter(unsigned int partDim, TunType type, bool isControl) 
  : SplitBase(partDim, type, isControl)
{
  m_addDir = false;
}
#endif

SplitSingleIter::~SplitSingleIter()
{
#if TWOD
  if (m_msizes) {
    delete [] m_msizes;
    m_msizes = NULL;
    delete [] m_nsizes;
    m_nsizes = NULL;
#if DODM
    delete [] m_mlsizes;
    m_mlsizes = NULL;
    delete [] m_nlsizes;
    m_nlsizes = NULL;
#endif
  }
#else
  if (m_sizes) {
    delete [] m_sizes;
    m_sizes = NULL;
    delete [] m_lsizes;
    m_lsizes = NULL;
  }
#endif
}

Node* SplitSingleIter::BlankInst() 
{ 
#if TWOD
  return new SplitSingleIter(LASTPARTDIR,LASTTUNNEL,false);
#else
  return new SplitSingleIter(99,LASTTUNNEL,false);
#endif
}

Name SplitSingleIter::GetName(ConnNum num) const 
{
  return GetName(num, GetLoopType());
}

Name SplitSingleIter::GetName(ConnNum num, LoopType type) const 
{
  Name name;
  if (type == BLISLOOP) {
    if (m_tunType == SETTUNOUT)
      return ((LoopTunnel*)Input(0))->GetOrigName();
    else if (m_tunType == POSSTUNOUT)
      return ((LoopTunnel*)Input(m_inputs.size()-1))->GetOrigName();
    else if (m_tunType == SETTUNIN) {
      name.m_name = GetInputName(0).str();
      if (name.m_name[name.m_name.length()-1] != '_')
        name.m_name = name.m_name + "_";
#if TWOD
  if (m_addDir) {
    switch(m_dir)
      {
      case (PARTDOWN):
	name.m_name += "PD_";
	break;
      case (PARTRIGHT):
	name.m_name += "PR_";
	break;
      case (PARTUPWARD):
	name.m_name += "PU_";
	break;
      case (PARTLEFT):     
	name.m_name += "PL_";
	break;
      case (PARTDIAG):
	name.m_name += "PDi_";
	break;
      case (PARTDIAGBACK):
	name.m_name += "PDB_";
	break;
      default:
	LOG_FAIL("replacement for throw call");
	throw;
      }
  }
#else
  if (m_addDir) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
#endif
      return name;
    }
  }
  name = GetInputName(0);

  if (type == BLISLOOP && name.m_name[name.m_name.length()-1] != '_') {
    //while constructing a new set, you put down a POSSTUNIN.
    // its input might be a POSSTUNIN (since the SETTUNIN is inserted
    // during set contruction). 
    // in such a case, the "_" would not have been added
    name.m_name = name.m_name + "_";
  }

	  

  if (m_tunType == POSSTUNIN) {
#if TWOD
    switch(m_dir)
      {
      case (PARTDOWN):
      case (PARTRIGHT):
      case (PARTUPWARD):
      case (PARTLEFT):     
        if (num == 0) 
          name.m_name += "0";
        else if (num == 1)
          name.m_name += "1";
        else if (num == 2)
          name.m_name += "2";
        else if (num == 3)
          break;
        else {
          LOG_FAIL("replacement for throw call");
	  throw;
	}
	break;
      case (PARTDIAG):
      case (PARTDIAGBACK):
        if (num == 0)
          name.m_name += "00";
        else if (num == 1)
          name.m_name += "10";
        else if (num == 2)
          name.m_name += "20";
        else if (num == 3)
          name.m_name += "01";
        else if (num == 4)
          name.m_name += "11";
        else if (num == 5)
          name.m_name += "21";
        else if (num == 6)
          name.m_name += "02";
        else if (num == 7)
          name.m_name += "12";
        else if (num == 8)
          name.m_name += "22";
        else if (num == 9)
          break;
        else {
          LOG_FAIL("replacement for throw call");
	  throw;
	}
	break;
      default:
        cout << "bad dir\n";
        LOG_FAIL("replacement for throw call");
	throw;
        break;
      }

    RealLoop *loop = (RealLoop*)(((Tunnel*)Input(0))->m_pset);
    if (loop->IsUnrolled())
      name.m_name += "_iter" + std::to_string((long long int) loop->GetCurrIter());
    return name;
#else
    if (m_partDim >= InputNumDims(0)) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    name.m_name += LoopLevel();
    if (num <= 2) 
      name.m_name += "_part" + std::to_string(m_partDim) + "_" + std::to_string(num);
    else if (num > 3) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
#endif
    return name;
  }
  else
    return name;
}

void SplitSingleIter::Prop()
{
  if (!IsValidCost(m_cost)) {
    LoopTunnel::Prop();


    if (m_tunType == POSSTUNIN 
	&& !IsSplit()) 
      {
	Node *found = NULL;
	NodeConnVecIter iter = m_children.begin();
	for( ; iter != m_children.end(); ++iter) {
	  if ((*iter)->m_num == (NumOutputs()-1)) {
	    if (found) {
	      LOG_FAIL("replacement for throw call");
	      throw;
	    } else {
	      found = (*iter)->m_n;
	    }
	  }
	}
      }

#if TWOD
    if ((m_dir == PARTDOWN || m_dir == PARTUPWARD)
	&& (GetUpStat(TL) != GetUpStat(TR) || GetUpStat(BL) != GetUpStat(BR)))
      {
	cout << "bad statuses\n";
	LOG_FAIL("replacement for throw call");
	throw;
      }
    else if ((m_dir == PARTRIGHT ||m_dir == PARTLEFT)
	     && (GetUpStat(TL) != GetUpStat(BL) || GetUpStat(TR) != GetUpStat(BR)))
      {
	cout << "bad statuses\n";
	LOG_FAIL("replacement for throw call");
	throw;
      }
#endif  
    if (m_tunType == POSSTUNIN) {
      for(ConnNum i = 0; i < m_inputs.size(); ++i) {
	if (Input(i)->GetNodeClass() != SplitSingleIter::GetClass()) {
	  LOG_FAIL("replacement for throw call");
	  throw;
	}
      }
    }
    else if (m_tunType == SETTUNIN) {
      if (!m_pset->IsLoop()) {
	if (m_pset->IsReal()) {
	  cout << "m_pset is real\n";
	} else {
	  cout << "m_pset is shadow\n";
	}
	cout << m_pset << endl;
	LOG_FAIL("replacement for throw call");
	throw;
      }
      if (m_inputs.size() != 1) {
	cout << "split has wrong number of inputs\n";
	LOG_FAIL("replacement for throw call");
	throw;
      }
      for (unsigned int i = 0; i < m_children.size(); ++i) {
	if (Child(i)->GetNodeClass() != SplitSingleIter::GetClass()) {
	  LOG_FAIL("replacement for throw call");
	  throw;
	}
      }
    }
    else {
      cout << "bad tunType\n";
      LOG_FAIL("replacement for throw call");
      throw;
    }

    m_cost = ZERO;
    if ((dynamic_cast<const LoopInterface*>(GetMyLoop()))->GetBS() == ZERO) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    if (m_tunType == POSSTUNIN) {
      if (Input(0)->GetNodeClass() != SplitSingleIter::GetClass()) {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      Node *child = NULL;
      unsigned int size = m_children.size();  
      for (unsigned int i = 0; i < size; ++i) {
#if TWOD
        if (m_children[i]->m_num == GetNumElems(m_dir)) {
#else
        if (m_children[i]->m_num == 3) {
#endif
          if (child) {
	    cout << child << endl;
	    cout << child->GetType() << endl;
	    cout << m_children[i]->m_n->GetType() << endl;
	    cout << m_children[i]->m_n << endl;
            LOG_FAIL("replacement for throw call");
	    throw;
	  }
          else
            child = m_children[i]->m_n;
        }
      }
      if (!child) {
        cout << GetNameStr(0) << endl;
        //m_poss->ForcePrint();
        LOG_FAIL("replacement for throw call");
	throw;
      }
      else if (child->GetNodeClass() != CombineSingleIter::GetClass()) {
        LOG_FAIL("replacement for throw call");
	throw;
      } else {
        CombineSingleIter* com = (CombineSingleIter*)child;
#if TWOD
        if (com->m_dir != m_dir) {
#else
	  if (com->m_partDim != m_partDim) {
#endif
          LOG_FAIL("replacement for throw call");
	  throw;
	  }
      }
    }
  }
}

#if TWOD
const Sizes* SplitSingleIter::GetM(ConnNum num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0) {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      return GetInputM(0);
    case (POSSTUNIN):
      if (num < GetNumElems(m_dir)) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_msizes) {
	  cout << Input(0)->GetType() << endl;
	  cout << GetMyLoop()->GetControl()->m_msizes << endl;
	  cout << GetMyLoop()->GetControl()->m_nsizes << endl;
	  cout << "loop " << m_pset << endl;
	  cout << "loop " << input->m_pset << endl;
	  
          LOG_FAIL("replacement for throw call");
	  throw;
	}
        return &(input->m_msizes[num]);
      }
      else if (num == GetNumElems(m_dir)) {
        return GetInputM(0);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

const Sizes* SplitSingleIter::GetN(ConnNum num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0) {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      return GetInputN(0);
    case (POSSTUNIN):
      if (num < GetNumElems(m_dir)) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_nsizes) {
	  cout << "Error: m_nsizes is null for POSSTUNIN" << endl;
          LOG_FAIL("replacement for throw call");
	  throw;
	}
        return &(input->m_nsizes[num]);
      }
      else if (num == GetNumElems(m_dir)) {
        return GetInputN(0);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

#if DODM
const Sizes* SplitSingleIter::LocalM(ConnNum num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0) {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      return InputLocalM(0);
    case (POSSTUNIN):
      if (num < GetNumElems(m_dir)) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_mlsizes) {
          LOG_FAIL("replacement for throw call");
	  throw;
	}
        return input->m_mlsizes+num;
      }
      else if (num == GetNumElems(m_dir)) {
        return InputLocalM(0);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

const Sizes* SplitSingleIter::LocalN(ConnNum num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0) {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      return InputLocalN(0);
    case (POSSTUNIN):
      if (num < GetNumElems(m_dir)) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_nlsizes) {
          LOG_FAIL("replacement for throw call");
	  throw;
	}
        return input->m_nlsizes+num;
      }
      else if (num == GetNumElems(m_dir)) {
        return InputLocalN(0);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}
#endif


void SplitSingleIter::GetSizes(ConnNum num, unsigned int numIters,
		     Size bs, unsigned int parFactor,
		     Size m, Size n,
		     Sizes &ms, Sizes &ns)
{
  if (m_tunType != SETTUNIN) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  //  if (GetLoopType() != ELEMLOOP) {
  //    cout << "check that the m and n values are actually local sizes\n";
  //    cout << "also use currlm = currm...\n";
  //    LOG_FAIL("replacement for throw call");
  //  }
  switch(m_dir)
    {
    case (PARTDOWN):
      if (num == 0) {
	ms.AddSizesWithLimit(0,bs,m, parFactor);
	ns.AddRepeatedSizes(n,numIters, parFactor);
      }
      else if (num ==1) {
	if (ceil(m/bs) != numIters) {
	  LOG_FAIL("replacement for throw call");
	  throw;
	}
	ms.AddMidSizes(bs, m, parFactor);
	ns.AddRepeatedSizes(n,numIters, parFactor);
      }
      else if (num == 2) {
	ms.AddSizesWithLimit(m,-1*bs,0, parFactor);
	ns.AddRepeatedSizes(n,numIters, parFactor);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      break;
    case (PARTUPWARD):
      if (num == 0) {
	ms.AddSizesWithLimit(m, -1*bs, 0, parFactor);
	ns.AddRepeatedSizes(n,numIters, parFactor);
      }
      else if (num ==1) {
	ms.AddMidSizes(bs,m, parFactor);
	ns.AddRepeatedSizes(n,numIters, parFactor);
      }
      else if (num == 2) {
	ms.AddSizesWithLimit(0,bs,m, parFactor);
	ns.AddRepeatedSizes(n,numIters, parFactor);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      break;
    case (PARTRIGHT):
      if (num == 0) {
	ms.AddRepeatedSizes(m,numIters, parFactor);
	ns.AddSizesWithLimit(0,bs,n, parFactor);
      }
      else if (num == 1) {
	ms.AddRepeatedSizes(m,numIters, parFactor);
	ns.AddMidSizes(bs,n, parFactor);
      }
      else if (num == 2) {
	ms.AddRepeatedSizes(m,numIters, parFactor);
	ns.AddSizesWithLimit(n,-1*bs,0, parFactor);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      break;
    case (PARTLEFT):
      if (num == 0) {
	ms.AddRepeatedSizes(m,numIters, parFactor);
	ns.AddSizesWithLimit(n,-1*bs,0, parFactor);
      }
      else if (num == 1) {
	ms.AddRepeatedSizes(m,numIters, parFactor);
	ns.AddMidSizes(bs,n, parFactor);
      }
      else if (num == 2) {
	ms.AddRepeatedSizes(m,numIters, parFactor);
	ns.AddSizesWithLimit(0,bs,n, parFactor);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      break;
    case (PARTDIAG):
      //First column
      if (num == 0) {
	ms.AddSizesWithLimit(0,bs,m, parFactor);
	ns.AddSizesWithLimit(0,bs,n, parFactor);
      }
      else if (num == 1) {
	ms.AddMidSizes(bs,m, parFactor);
	ns.AddSizesWithLimit(0,bs,n, parFactor);
      }
      else if (num == 2) {
	ms.AddSizesWithLimit(m,-1*bs,0, parFactor);
	ns.AddSizesWithLimit(0,bs,n, parFactor);
      }
      //Second column
      else if (num == 3) {
	ms.AddSizesWithLimit(0,bs,m, parFactor);
	ns.AddMidSizes(bs,n, parFactor);
      }
      else if (num == 4) {
	ms.AddMidSizes(bs,m, parFactor);
	ns.AddMidSizes(bs,n, parFactor);
      }
      else if (num == 5) {
	ms.AddSizesWithLimit(m,-1*bs,0, parFactor);
	ns.AddMidSizes(bs,n, parFactor);
      }
      //Third column
      else if (num == 6) {
	ms.AddSizesWithLimit(0,bs,m, parFactor);
	ns.AddSizesWithLimit(n,-1*bs,0, parFactor);
      }
      else if (num == 7) {
	ms.AddMidSizes(bs,m, parFactor);
	ns.AddSizesWithLimit(n,-1*bs,0, parFactor);
      }
      else if (num == 8) {
	ms.AddSizesWithLimit(m,-1*bs,0, parFactor);
	ns.AddSizesWithLimit(n,-1*bs,0, parFactor);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      break;
    case (PARTDIAGBACK):
      //First column
      if (num == 0) {
	ms.AddSizesWithLimit(m,-1*bs,0, parFactor);
	ns.AddSizesWithLimit(n,-1*bs,0, parFactor);
      }
      else if (num == 1) {
	ms.AddMidSizes(bs,m, parFactor);
	ns.AddSizesWithLimit(n,-1*bs,0, parFactor);
      }
      else if (num == 2) {
	ms.AddSizesWithLimit(0,bs,m, parFactor);
	ns.AddSizesWithLimit(n,-1*bs,0, parFactor);
      }
      //Second column
      else if (num == 3) {
	ms.AddSizesWithLimit(m,-1*bs,0, parFactor);
	ns.AddMidSizes(bs,n, parFactor);
      }
      else if (num == 4) {
	ms.AddMidSizes(bs,m, parFactor);
	ns.AddMidSizes(bs,n, parFactor);
      }
      else if (num == 5) {
	ms.AddSizesWithLimit(0,bs,m, parFactor);
	ns.AddMidSizes(bs,n, parFactor);
      }
      //Third column
      else if (num == 6) {
	ms.AddSizesWithLimit(m,-1*bs,0, parFactor);
	ns.AddSizesWithLimit(0,bs,n, parFactor);
      }
      else if (num == 7) {
	ms.AddMidSizes(bs,m, parFactor);
	ns.AddSizesWithLimit(0,bs,n, parFactor);
      }
      else if (num == 8) {
	ms.AddSizesWithLimit(0,bs,m, parFactor);
	ns.AddSizesWithLimit(0,bs,n, parFactor);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      break;
    default:
      LOG_FAIL("replacement for throw call");
      throw;
      break;
    }
}
#else
 const Dim SplitSingleIter::NumDims(ConnNum num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0) {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      return InputNumDims(0);
    case (POSSTUNIN):
      if (num > 3) {
	LOG_FAIL("replacement for throw call");
	throw;
      }
      return InputNumDims(0);
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

const Sizes* SplitSingleIter::Len(ConnNum num, Dim dim) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0) {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      return InputLen(0,dim);
    case (POSSTUNIN):
      if (num < 3) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_sizes) {
          LOG_FAIL("replacement for throw call");
	  throw;
	}
	if (dim < m_partDim)
	  return &(input->m_sizes[dim]);
	else if (dim == m_partDim)
	  return &(input->m_sizes[dim+num]);
	else
	  return &(input->m_sizes[dim+2]);
      }
      else if (num == 3) {
        return InputLen(0,dim);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

const Sizes* SplitSingleIter::LocalLen(ConnNum num, Dim dim) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0) {
        LOG_FAIL("replacement for throw call");
	throw;
      }
      return InputLocalLen(0,dim);
    case (POSSTUNIN):
      if (num < 3) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_lsizes) {
          LOG_FAIL("replacement for throw call");
	  throw;
	}
	if (dim < m_partDim)
	  return &(input->m_lsizes[dim]);
	else if (dim == m_partDim)
	  return &(input->m_lsizes[dim+num]);
	else
	  return &(input->m_lsizes[dim+2]);
      }
      else if (num == 3) {
        return InputLocalLen(0,dim);
      }
      else {
        LOG_FAIL("replacement for throw call");
	throw;
      }
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}
#endif

#if TWOD
Tunnel* SplitSingleIter::GetSetTunnel()
{
  SplitSingleIter *tun;
  if (m_tunType == POSSTUNIN || m_tunType == SETTUNIN)
    tun = new SplitSingleIter(m_dir, SETTUNIN, m_isControlTun);
  else if (m_tunType == POSSTUNOUT || m_tunType == SETTUNOUT)
    tun = new SplitSingleIter(m_dir, SETTUNOUT, m_isControlTun);
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  tun->CopyTunnelInfo(this);
  return tun;
}
#else
Tunnel* SplitSingleIter::GetSetTunnel()
{
  SplitSingleIter *tun;
  if (m_tunType == POSSTUNIN)
    tun = new SplitSingleIter(m_partDim, SETTUNIN, m_isControlTun);
  else if (m_tunType == POSSTUNOUT)
    tun = new SplitSingleIter(m_partDim, SETTUNOUT, m_isControlTun);
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  tun->CopyTunnelInfo(this);
  return tun;
}
#endif

void SplitSingleIter::PrintCode(IndStream &out)
{
  if (m_tunType != POSSTUNIN)
    return;
#if TWOD
  LoopType loopType = GetLoopType();
  if (loopType == ELEMLOOP) {
    switch(m_dir) {
    case (PARTDOWN):
      out.Indent();
      *out << "RepartitionDown\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(0) << "T,  " << GetNameStr(0) << ",\n"
	   << out.Tabs(0)
	   << "  /**/ /**/\n"
	   << out.Tabs(0)
	   << "       " << GetNameStr(1) << ",\n"
	   << out.Tabs(0)
	   << "  " << GetInputNameStr(0) << "B, " << GetNameStr(2) << " );\n";
      break;
    case (PARTUPWARD):
      out.Indent();
      *out << "RepartitionUp\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(0) << "T,  " << GetNameStr(0) << ",\n"
	   << out.Tabs(0)
	   << "       " << GetNameStr(1) << ",\n"
	   << out.Tabs(0)
	   << "  /**/ /**/\n"
	   << out.Tabs(0)
	   << "  " << GetInputNameStr(0) << "B, " << GetNameStr(2) << " );\n";
      break;
    case (PARTRIGHT):
      out.Indent();
      *out << "RepartitionRight\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(0) << "L, /**/     " << GetInputNameStr(0) << "R,\n"
	   << out.Tabs(0)
	   << "( " << GetNameStr(0) << ", /**/ " << GetNameStr(1) << ", " << GetNameStr(2) << " );\n";
      break;
    case (PARTLEFT):
      out.Indent();
      *out << "RepartitionLeft\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(0) << "L, /**/     " << GetInputNameStr(0) << "R,\n"
	   << out.Tabs(0)
	   << "( " << GetNameStr(0) << ", " << GetNameStr(1) << ", /**/ " << GetNameStr(2) << " );\n";
      break;
    case (PARTDIAG):
      out.Indent();
      *out << "RepartitionDownDiagonal\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(0) << "TL, /**/ " << GetInputNameStr(0) << "TR,  "
	   << GetNameStr(0) << ", /**/ " << GetNameStr(3) << ", " << GetNameStr(6) << ",\n"
	   << out.Tabs(0)
	   << " /*************/ /******************/\n"
	   << out.Tabs(0)
	   << "       /**/       "
	   << GetNameStr(1) << ", /**/ " << GetNameStr(4) << ", " << GetNameStr(7) << ",\n"
	   << out.Tabs(0)
	   << "  " << GetInputNameStr(0) << "BL, /**/ " << GetInputNameStr(0) << "BR,  "
	   << GetNameStr(2) << ", /**/ " << GetNameStr(5) << ", " << GetNameStr(8) << " );\n";
      break;
    case (PARTDIAGBACK):
      out.Indent();
      *out << "RepartitionLeftDiagonal\n"
	   << out.Tabs(0)
	   << "( " << GetInputNameStr(0) << "TL, /**/ " << GetInputNameStr(0) << "TR,  "
	   << GetNameStr(0) << ", " << GetNameStr(3) << ", /**/ " << GetNameStr(6) << ",\n"
	   << out.Tabs(0)
	   << "       /**/       "
	   << GetNameStr(1) << ", " << GetNameStr(4) << ", /**/ " << GetNameStr(7) << ",\n"
	   << out.Tabs(0)
	   << " /*************/ /******************/\n"
	   << out.Tabs(0)
	   << "  " << GetInputNameStr(0) << "BL, /**/ " << GetInputNameStr(0) << "BR,  "
	   << GetNameStr(2) << ", " << GetNameStr(5) << ", /**/ " << GetNameStr(8) << " );\n";
      break;
    default:
      cout << "bad dir\n";
      break;
    }
  }
  else if (loopType == BLISLOOP) {
    Node *setTunIn = Input(0);
    string inName = setTunIn->Input(0)->GetNameStr(setTunIn->InputConnNum(0));
    for (unsigned int num = 0; num < NumOutputs(); ++num) {
      bool hasPrinted = false;
      NodeConnVecIter iter = m_children.begin();
      for (; !hasPrinted && iter != m_children.end(); ++iter) {
	if ((*iter)->m_num != num)
	  continue;
	Node *child = (*iter)->m_n;
	if (!child->IsTunnel(POSSTUNOUT)) {
	  hasPrinted = true;
	  out.Indent();
	  *out << "obj_t " << GetNameStr(num) << ";\n";	
	  out.Indent();
	  switch(m_dir) {
	  case (PARTUPWARD):
	    *out << "bli_acquire_mpart_b2t";
	    break;
	  case (PARTDOWN):
	    *out << "bli_acquire_mpart_t2b";
	    break;
	  case (PARTLEFT):
	    *out << "bli_acquire_mpart_r2l";
	    break;
	  case (PARTRIGHT):
	    *out << "bli_acquire_mpart_l2r";
	    break;
	  case (PARTDIAG):
	    *out << "bli_acquire_mpart_tl2br";
	    break;
	  case (PARTDIAGBACK):
	    *out << "bli_acquire_mpart_br2tl";
	    break;
	  default:
	    LOG_FAIL("replacement for throw call");
	    throw;
	  }
	  string part;
	  switch(m_dir) {
	  case (PARTUPWARD):
	  case (PARTDOWN):
	  case (PARTLEFT):
	  case (PARTRIGHT):
	    if (num == 0)
	      part = "BLIS_SUBPART0";
	    else if (num == 1)
	      part = "BLIS_SUBPART1";
	    else if (num == 2)
	      part = "BLIS_SUBPART2";
	    else {
	      LOG_FAIL("replacement for throw call");
	      throw;
	    }
	    break;
	  case (PARTDIAG):
	  case (PARTDIAGBACK):
	    if (num == 0)
	      part = "BLIS_SUBPART00";
	    else if (num == 1)
	      part = "BLIS_SUBPART10";
	    else if (num == 2)
	      part = "BLIS_SUBPART20";
	    else if (num == 3)
	      part = "BLIS_SUBPART01";
	    else if (num == 4)
	      part = "BLIS_SUBPART11";
	    else if (num == 5)
	      part = "BLIS_SUBPART21";
	    else if (num == 6)
	      part = "BLIS_SUBPART02";
	    else if (num == 7)
	      part = "BLIS_SUBPART12";
	    else if (num == 8)
	      part = "BLIS_SUBPART22";
	    else {
	      LOG_FAIL("replacement for throw call");
	      throw;
	    }
	    break;
	  default:
	    LOG_FAIL("replacement for throw call");
	    throw;
	  }
	  string loopLevel = out.LoopLevel();
	  *out << "( " << part << ", idx" << loopLevel << ", bs" << loopLevel
	       << ", &" << inName << ", &" << GetNameStr(num) << " );\n";
	}
      }
    }
  }
#elif DOTENSORS
  Name nameT = GetInputName(0);
  nameT.m_name += LoopLevel();
  Name nameB = nameT;
  nameT.m_name += "_part" + std::to_string(m_partDim) + "T";
  nameB.m_name += "_part" + std::to_string(m_partDim) + "B";
    
  out.Indent();
  *out << "RepartitionDown\n"
       << out.Tabs(0)
       << "( " << nameT.str() << ",  " << GetNameStr(0) << ",\n"
       << out.Tabs(0)
       << "  /**/ /**/\n"
       << out.Tabs(0)
       << "       " << GetNameStr(1) << ",\n"
       << out.Tabs(0)
       << "  " << nameB.str() << ", " << GetNameStr(2) << ", " 
       <<  m_partDim << ", " 
       << GetMyLoop()->GetBS() << " );\n";

#else
  *out << "need split print code\n";
  LOG_FAIL("replacement for throw call");
  throw;
#endif
}

void SplitSingleIter::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  SplitBase::Duplicate(orig, shallow, possMerging);
  const SplitSingleIter *split = (SplitSingleIter*)orig;
  m_addDir = split->m_addDir;
#if DOLLDLA
  m_info = split->m_info;
#endif
}

NodeType SplitSingleIter::GetType() const
{
#if TWOD
  return "SplitSingleIter " + PartDirToStr(m_dir)  + "( " + Tunnel::GetType() + " )";
#else
  string tmp = "SplitSingleIter";
  tmp += m_partDim;
  tmp += "_" + LoopLevel();
  return tmp  + "( " + Tunnel::GetType() + " )";
#endif
}

unsigned int SplitSingleIter::NumOutputs() const 
{
#if TWOD
  return GetNumElems(m_dir)+1;
#else
  return 4;
#endif
}

bool SplitSingleIter::QuadInUse(Quad quad, bool atEnd) const
{
  if (m_pset && !m_pset->IsReal())
    return ((SplitSingleIter*)GetRealTunnel())->QuadInUse(quad,atEnd);
  if (m_tunType == SETTUNIN) {
    Node *child = Child(0);
    if (!child->IsLoopTunnel()) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    return ((LoopTunnel*)child)->QuadInUse(quad, atEnd);
  }
  else if (m_tunType == POSSTUNIN) {
    NodeConnVecConstIter iter = m_children.begin();
    for( ; iter != m_children.end(); ++iter) {
      bool check = false;
      ConnNum num = (*iter)->m_num;
#if DOTENSORS
      switch (PARTDOWN) {
#else
      switch (m_dir) {
#endif
      case (PARTDOWN):
	if ((quad == TL || quad == TR) && ((!atEnd && num == 0) || (atEnd && (num == 0 || num == 1))))
	  check = true;
	else if ((quad == BL || quad == BR) && ((!atEnd && (num == 1 || num == 2)) || (atEnd && (num == 2))))
	  check = true;
	else
	  check = false;
	break;
      case (PARTUPWARD):
	if ((quad == TL || quad == TR) && ((!atEnd && (num == 0 || num == 1)) || (atEnd && (num == 0))))
	  check = true;
	else if ((quad == BL || quad == BR) && ((!atEnd && num == 2) || (atEnd && (num == 1 || num == 2))))
	  check = true;
	else
	  check = false;
	break;
      case (PARTRIGHT):
	if ((quad == TL || quad == BL) && ((!atEnd && num == 0) || (atEnd && (num == 0 || num == 1))))
	  check = true;
	else if ((quad == TR || quad == BR) && ((!atEnd && (num == 1 || num == 2)) || (atEnd && num == 2)))
	  check = true;
	else
	  check = false;
	break;
      case (PARTLEFT):
	if ((quad == TL || quad == BL) && ((!atEnd && (num == 0 || num == 1)) || (atEnd && num == 0)))
	  check = true;
	else if ((quad == TR || quad == BR) && ((!atEnd && num == 2) || (atEnd && (num == 1 || num == 2))))
	  check = true;
	else
	  check = false;
	break;
      case (PARTDIAG):
	if (quad == TL && ((!atEnd && num == 0) || (atEnd && (num == 0 || num == 1 || num == 3 || num == 4))))
	  check = true;
	else if (quad == BL && ((!atEnd && (num == 1 || num == 2)) || (atEnd && (num == 2 || num == 5))))
	  check = true;
	else if (quad == TR && ((!atEnd && (num == 3 || num == 6)) || (atEnd && (num == 6 || num == 7))))
	  check = true;
	else if (quad == BR && ((!atEnd && (num == 4 || num == 5 || num == 7 || num == 8)) || (atEnd && num == 8)))
	  check = true;
	else
	  check = false;
	break;
      case (PARTDIAGBACK):
	if (quad == TL && ((!atEnd && (num == 0 || num == 1 || num == 3 || num == 4)) || (atEnd && num == 0)))
	  check = true;
	else if (quad == BL && ((!atEnd && (num == 2 || num == 5)) || (atEnd && (num == 1 || num == 2))))
	  check = true;
	else if (quad == TR && ((!atEnd && (num == 6 || num == 7)) || (atEnd && (num == 3 || num == 6))))
	  check = true;
	else if (quad == BR && ((!atEnd && (num == 8)) || (atEnd && (num == 4 || num == 5 || num == 7 || num == 8))))
	  check = true;
	else
	  check = false;
	break;
      default:
	LOG_FAIL("replacement for throw call");
	throw;
      }
      if (check && !(*iter)->m_n->IsTunnel(POSSTUNOUT))
        return true;
    }
    return false;
  }
    else {
    LOG_FAIL("replacement for throw call");
    throw;
    }
}

void SplitSingleIter::PrintVarDeclarations(BSSize bs, IndStream &out) const
{
  //update PrintIncrementAtEndOfLoop, too
#if DOLLDLA
  if (m_tunType != POSSTUNIN) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  const string name = GetInputNameStr(0);
  if (m_dir != PARTDOWN && m_dir != PARTRIGHT) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (PartInUse(0) || PartInUse(2)) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  out.Indent();

  BasePSet *pset = ((Tunnel*)Input(0))->m_pset;
  if (pset->IsReal() && ((RealLoop*)pset)->IsUnrolled()) {
    unsigned int iter = ((RealLoop*)pset)->GetCurrIter();
    *out << LLDLAPartVarName(name, 1) << "_iter"
	 << iter
	 << " = " << name;
    
    const DataTypeInfo &type = InputDataType(0);
    *out << " + " << iter << " * (";
    if (m_dir == PARTDOWN) {
      if (!IsUnitStride(type.m_rowStride))
	*out << type.m_rowStrideVar << " * ";
      *out << bs.VarName();
    }
    else if (m_dir == PARTRIGHT) {
      if (!IsUnitStride(type.m_colStride))
	*out << type.m_colStrideVar << " * ";
      *out << bs.VarName();
    }
    else {
      LOG_FAIL("replacement for throw call");
      throw;
    }

    *out << ");\n";
  }
  else {
    *out << LLDLAPartVarName(name, 1)
	 << " = " << name << ";\n";
  }
#elif DOTENSORS
  Name nameT = GetInputName(0);
  nameT.m_name += LoopLevel();
  Name nameB = nameT;
  nameT.m_name += "_part" + std::to_string(m_partDim) + "T";
  nameB.m_name += "_part" + std::to_string(m_partDim) + "B";
  out.Indent();
  *out << "PartitionDown(" << GetInputNameStr(0) << ", " 
       << nameT.str() << ", "
       << nameB.str() << ", "
       << m_partDim << ", 0);\n";
#endif
}

bool SplitSingleIter::PartInUse(unsigned int partNum) const
{
  if (m_tunType != POSSTUNIN) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  NodeConnVecConstIter iter = m_children.begin();
  for(; iter != m_children.end(); ++iter) {
    const NodeConn *conn = *iter;
    if (conn->m_num == partNum) {
      if (!conn->m_n->IsTunnel(POSSTUNOUT))
	return true;
    }
  }
  return false;	    
}

void SplitSingleIter::AddVariables(VarSet &set) const
{
#if DOLLDLA
  if (m_tunType != POSSTUNIN)
    return;
  if (GetLoopType() != LLDLALOOP) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (m_dir != PARTDOWN && m_dir != PARTRIGHT) {
    LOG_FAIL("replacement for throw call");
    throw;
  }

  const string name = GetInputNameStr(0);

  BasePSet *pset = ((Tunnel*)Input(0))->m_pset;
  if (GetDataType() == REAL_SINGLE) {
    //    cout << "Var " << name << " has type is float\n";
  } else if (GetDataType() == REAL_DOUBLE) {
    //    cout << "Var " << name << " has type is double\n";
  } else {
    cout << "Error: Bad type, var " << name << " in SplitSingleIter::AddVariables\n";
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (!pset->IsReal() || !((RealLoop*)pset)->IsUnrolled()) {
    if (PartInUse(0)) {
      Var var(name, 0, GetDataType());
      set.insert(var);
    }
    Var var1(name, 1, GetDataType());
    set.insert(var1);
    if (PartInUse(2)) {
      Var var(name, 2, GetDataType());
      set.insert(var);
    }

    if (m_isControlTun) {
      string loopLevel = GetLoopLevel(-1);
      Var var(DirectVarDeclType, "int lcv"+loopLevel+";\n", GetDataType());
      set.insert(var);
    }
  }
  else if (pset->IsReal()) {
    if (PartInUse(0) || PartInUse(2)) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    unsigned int numIters = NumIters(0);
    for(unsigned int i = 0; i < numIters; ++i) {
      if (GetDataType() == REAL_DOUBLE) {
	Var varD(DirectVarDeclType, "double *" + LLDLAPartVarName(name, 1) + "_iter" + std::to_string((long long int) i) + ";", GetDataType());
	set.insert(varD);
      } else if (GetDataType() == REAL_SINGLE) {
	Var varS(DirectVarDeclType, "float *" + LLDLAPartVarName(name, 1) + "_iter" + std::to_string((long long int) i) + ";", GetDataType());
	set.insert(varS);
      }
      }
      }
#elif DOTENSORS
  {
    Name name = GetInputName(0);
    name.m_name += LoopLevel();
    name.m_name += "_part" + std::to_string(m_partDim) + "T";
    Var var(name);
    set.insert(var);
  }
{
    Name name = GetInputName(0);
    name.m_name += LoopLevel();
    name.m_name += "_part" + std::to_string(m_partDim) + "B";
    Var var(name);
    set.insert(var);
  }
  {
    Var var(GetName(0));
    set.insert(var);
  }
  {
    Var var(GetName(1));
    set.insert(var);
  }
  {
    Var var(GetName(2));
    set.insert(var);
  }
#endif
}

CombineSingleIter* SplitSingleIter::CreateMatchingCombine(int numArgs, ...)
{
  int numComIns = NumOutputs();
  int j = 0;
#if TWOD
  CombineSingleIter *com = new CombineSingleIter(m_dir, POSSTUNOUT);
#else
  CombineSingleIter *com = new CombineSingleIter(m_partDim, POSSTUNOUT);
#endif
  
  va_list listPointer;
  va_start (listPointer, numArgs);
  for(int i = 0; i < numArgs; i++) {
    int inputNum = va_arg(listPointer, int);
    Node *node = va_arg(listPointer, Node* );
    int num = va_arg(listPointer, int);
    while (j < inputNum && j < numComIns) {
      com->AddInput(this, j);
      ++j;
    }
    if (j >= numComIns) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    com->AddInput(node, num);
    ++j;
  }
  while (j < numComIns) {
    com->AddInput(this, j);
    ++j;
  }
  com->CopyTunnelInfo(this);
  return com;
}

#if TWOD
unsigned int SplitSingleIter::NumIters(Size bs, Size m, Size n) const
{
  if (!bs) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  switch (m_dir) {
    case(PARTDOWN):
    case(PARTDIAG):
    case(PARTUPWARD):
    case(PARTDIAGBACK):
      return ceil(((double)m)/bs);
    case(PARTLEFT):
    case(PARTRIGHT):
      return ceil(((double)n)/bs);
    case(LASTPARTDIR):
      LOG_FAIL("replacement for throw call");
      throw;
  }
  LOG_FAIL("replacement for throw call");
  throw;
}
#else
unsigned int SplitSingleIter::NumIters(Size bs, Size size) const
{
  if (!bs) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  return ceil(((double)size)/bs);
}
#endif

#if TWOD
unsigned int SplitSingleIter::NumIters(unsigned int iterNum) const
{
  Size bs = GetMyLoop()->GetBS();

  const Sizes *ms = GetInputM(0);
  const Sizes *ns = GetInputN(0);
  if (!ms || !ns) {
    if (Input(0)->m_flags & NODEBUILDFLAG)
      cout << "has built\n";
    else
      cout << "hasn't built\n";
    cout << Input(0)->GetNodeClass() << endl;
    cout << Input(0) << endl;
    LOG_FAIL("replacement for throw call");
    throw;
  }
  const Size m = (*ms)[iterNum];
  const Size n = (*ns)[iterNum];
  return NumIters(bs, m, n);
}
#else
unsigned int SplitSingleIter::NumIters(unsigned int iterNum) const
{
  Size bs = GetMyLoop()->GetBS();
  const Sizes *sizes = InputLen(0,m_partDim);
  if (!sizes) {
    if (Input(0)->m_flags & NODEBUILDFLAG)
      cout << "has built\n";
    else
      cout << "hasn't built\n";
    cout << Input(0)->GetNodeClass() << endl;
    cout << Input(0) << endl;
    LOG_FAIL("replacement for throw call");
    throw;
  }
  const Size m = (*sizes)[iterNum];
  return NumIters(bs,m);
}
#endif

void SplitSingleIter::FlattenCore(ofstream &out) const
{
  SplitBase::FlattenCore(out);
  WRITE(m_addDir);
}


void SplitSingleIter::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  SplitBase::UnflattenCore(in,info);
  READ(m_addDir);
}


unsigned int SplitSingleIter::NumberOfLoopExecs() const
{
#if TWOD
  unsigned int one = GetInputM(0)->NumSizes();
  unsigned int two = GetInputN(0)->NumSizes();

  if (one != two)
    throw;
  return one;
#else
  return InputLen(0,0)->NumSizes();
#endif
}
#if TWOD
void SplitSingleIter::StartFillingSizes()
{
  if (m_msizes) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (m_tunType != SETTUNIN)
    return;
  if (m_pset && !m_pset->IsReal())
    return;
  unsigned int numElems = GetNumElems(m_dir);
  m_msizes = new Sizes[numElems];
  m_nsizes = new Sizes[numElems];
#if DODM
  m_mlsizes = new Sizes[numElems];
  m_nlsizes = new Sizes[numElems];
#endif
  /*
#if DOLLDLA
  if (m_tunType == SETTUNIN) {
    m_info = InputDataType(0);
    switch (m_dir) {
    case (PARTDOWN):
      m_info.m_numRowsVar = "numRows" + GetLoopLevel();
      break;
    case (PARTRIGHT):
      m_info.m_numColsVar = "numCols" + GetLoopLevel();
      break;
    default:
      LOG_FAIL("replacement for throw call");
    }
  }
#endif
  */
}
#else
void SplitSingleIter::StartFillingSizes()
{
  if (m_sizes) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (m_tunType != SETTUNIN)
    return;
  if (m_pset && !m_pset->IsReal())
    return;
  unsigned int numDims = InputNumDims(0);
  //Num dims of sizes, but the m_partDim dimension has
  // 3 outputs
  m_sizes = new Sizes[numDims+2];
  m_lsizes = new Sizes[numDims+2];
}
#endif

void SplitSingleIter::ClearDataTypeCache()
{
#if TWOD
  if (!m_msizes)
    return;
  delete [] m_msizes;
  m_msizes = NULL;
  delete [] m_nsizes;
  m_nsizes = NULL;
#if DODM
  delete [] m_mlsizes;
  m_mlsizes = NULL;
  delete [] m_nlsizes;
  m_nlsizes = NULL;
#endif
#else
  if (!m_sizes)
    return;
  delete  [] m_sizes;
  m_sizes = NULL;
  delete [] m_lsizes;
  m_lsizes = NULL;
#endif
}

#if TWOD&&DODM
void SplitSingleIter::UpdateLocalSizes()
{
  if (m_pset && !m_pset->IsReal())
    return;
  const unsigned int numElems = GetNumElems(m_dir);
  const LoopType loopType = GetMyLoop()->GetType();
  if (loopType == ELEMLOOP) {
    const DistType t = InputDataType(0).m_dist;

    for (unsigned int subMat = 0; subMat < numElems; ++subMat) {
      GetLocalSizes(t, m_msizes+subMat, m_nsizes+subMat, m_mlsizes[subMat], m_nlsizes[subMat]);
    }
  }
  else {
    cout << loopType << endl;
    LOG_FAIL("replacement for throw call");
    throw;
  }
}
#elif DOTENSORS
void SplitSingleIter::UpdateLocalSizes()
{
  if (m_pset && !m_pset->IsReal())
    return;
  Dim numDims = InputNumDims(0);
  const DistType t = InputDataType(0).GetEffectiveDist();
  for (Dim dim = 0; dim < numDims; ++ dim) {
    if (dim < m_partDim)
      GetLocalSizes(t, dim, m_sizes+dim, m_lsizes+dim);
    else if (dim == m_partDim) {
      GetLocalSizes(t, dim, m_sizes+dim, m_lsizes+dim);
      GetLocalSizes(t, dim, m_sizes+dim+1, m_lsizes+dim+1);
      GetLocalSizes(t, dim, m_sizes+dim+2, m_lsizes+dim+2);
    }
    else 
      GetLocalSizes(t, dim, m_sizes+dim+2, m_lsizes+dim+2);
  }
}
#endif

#if DOLLDLA
string SplitSingleIter::LoopBound()
{
  const DataTypeInfo &info = InputDataType(0);
  switch(m_dir)
    {
    case (PARTDOWN):
      {
	string tmp = "(" + GetInputNameStr(0) +
	  " + " + info.m_numRowsVar;
	if (!IsUnitStride(info.m_rowStride))
	  tmp += " * " + info.m_rowStrideVar;
	return tmp + ")";
	break;
      }
    case (PARTRIGHT):
      {
	string tmp = "(" + GetInputNameStr(0) +
	  " + " + info.m_numColsVar;
	if (!IsUnitStride(info.m_colStride))
	  tmp += " * " + info.m_colStrideVar;
	return tmp + ")";
	break;
      }
    case (PARTDIAG):
      LOG_FAIL("replacement for throw call");
      throw;
      break;
    case (PARTUPWARD):
      LOG_FAIL("replacement for throw call");
      throw;
      break;
    case (PARTLEFT):     
      LOG_FAIL("replacement for throw call");
      throw;
      break;

    case (PARTDIAGBACK):
      LOG_FAIL("replacement for throw call");
      throw;
      break;
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}
#endif


void SplitSingleIter::PrintIncrementAtEndOfLoop(BSSize bs, IndStream &out) const
{
  // Update PrintVarDeclarations, too
  if (m_tunType != POSSTUNIN) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
#if DOLLDLA
  out.Indent(1);
  const DataTypeInfo &type = InputDataType(0);
  *out << LLDLAPartVarName(GetInputNameStr(0),1) << " += ";
  if (m_dir == PARTDOWN) {
    if (!IsUnitStride(type.m_rowStride))
      *out << type.m_rowStrideVar << " * ";
    *out << bs.VarName() << ";\n";
  }
  else if (m_dir == PARTRIGHT) {
    if (!IsUnitStride(type.m_colStride))
      *out << type.m_colStrideVar << " * ";
    *out << bs.VarName() << ";\n";
  }
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (PartInUse(0) || PartInUse(2)) {
    LOG_FAIL("replacement for throw call");
    throw;
  }

#endif
}

#if DOLLDLA
void SplitSingleIter::BuildDataTypeCache()
{
  SplitBase::BuildDataTypeCache();
  if (m_pset && !m_pset->IsReal())
    return;
  if (m_tunType == SETTUNIN) {
    m_info = InputDataType(0);
    switch (m_dir) {
    case (PARTDOWN):
      m_info.m_numRowsVar = "numRows" + GetNameStr(1);
      break;
    case (PARTRIGHT):
      m_info.m_numColsVar = "numCols" + GetNameStr(1);
      break;
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
  }
}


const DataTypeInfo& SplitSingleIter::DataType(ConnNum num) const
{
  if (m_tunType == SETTUNIN) {
    if (m_pset && !m_pset->IsReal()) {
      return GetRealTunnel()->DataType(num);
    }
    else {
      unsigned int numElems = GetNumElems(m_dir);
      if (num < numElems) {
	return m_info;
      }
      else if (num == numElems) {
	return InputDataType(0);
      }
      else {
	LOG_FAIL("replacement for throw call");
	throw;
      }
    }
  }
  else {
    return Input(0)->DataType(num);
  }
}
#endif


#if DOLLDLA
void SplitSingleIter::MigrateFromOldTun(Tunnel *tun)
{
  LoopTunnel::MigrateFromOldTun(tun);
  if (!tun->IsSplit()) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  m_info = ((SplitSingleIter*)tun)->m_info;
}
#endif


LoopTunnel* SplitSingleIter::GetMatchingOutTun() const
{
  if (m_tunType == SETTUNIN)
    if (m_pset->IsReal()) 
      return (LoopTunnel*)(((LoopTunnel*)Child(0))->GetMatchingOutTun()->Child(0));
    else {
      LoopTunnel *possTunOut = ((LoopTunnel*)(GetRealTunnel()->Child(0)))->GetMatchingOutTun();
      LoopTunnel *realSetTunOut = (LoopTunnel*)(possTunOut->Child(0));
      return (LoopTunnel*)(m_pset->m_outTuns[FindInTunVec(m_pset->GetReal()->m_outTuns, realSetTunOut)]);
    }
  else if (m_tunType != POSSTUNIN) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  

  //BAMTODO : instead of doing this search, the child connected to the 
  // second (1th) output should be the poss out tun
  NodeConnVecConstIter iter = m_children.begin();
  for(; iter != m_children.end(); ++iter) {
    const NodeConn *con = *iter;
    const Node *child = con->m_n;
    if (child->IsTunnel(POSSTUNOUT) &&
        ((Tunnel*)child)->IsCombine()) {
      return (LoopTunnel*)child;
    }
  }
  cout << "Didn't find matching out tun\n";
  LOG_FAIL("replacement for throw call");
  throw;
}


//If we're using this to get alignment information
// during transformation application, this poss
// won't be part of the PSet yet, so we have to 
// use node connections to work up/out
string SplitSingleIter::LoopLevel() const
{
  Poss *poss = m_poss;
  int level = 0;
  while (poss) {
    if (poss->m_inTuns.empty())
      break;
    Node *in = poss->m_inTuns[0];
    if (!in)
      break;
    if (in->IsLoopTunnel())
      ++level;
    in = in->Input(0);
    if (!in) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    if (poss == in->m_poss) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    poss = in->m_poss;
  }
  return (string)"_lvl" + ((char)(level+48));
}

#if TWOD
void SplitSingleIter::BuildSizes(bool buildCache, vector<int> &numItersVec, unsigned int parFactor)
{
  if (m_tunType != SETTUNIN)
    return;
  if (m_pset && !m_pset->IsReal())
    return;
  if (!m_msizes) {
    LOG_FAIL("replacement for throw call");
    throw;
  }

  const Sizes *ms = GetInputM(0);
  const Sizes *ns = GetInputN(0);

  if (!ms || !ns) {
    LOG_FAIL("replacement for throw call");
    throw;
  }

  unsigned int length = ms->NumSizes();
  unsigned int length2 = ns->NumSizes();

  if (length != length2)
      LOG_FAIL("replacement for throw call");

  unsigned int numExecs;
  if (buildCache) {
    if (!m_isControlTun)
      LOG_FAIL("replacement for throw call");
    numExecs = length;
    numItersVec.reserve(numExecs);
  }
  else {
    if (m_isControlTun)
      LOG_FAIL("replacement for throw call");
    numExecs = numItersVec.size();
    if (numExecs != length)
      LOG_FAIL("replacement for throw call");
  }

  const Size bs = GetMyLoop()->GetBS();
  const unsigned int numElems = GetNumElems(m_dir);
  
  bool foundOne = false;
  for (unsigned int subMat = 0; subMat < numElems; ++subMat) {
    bool found = false;
    NodeConnVecIter tunIter = m_children.begin();
    for(; tunIter != m_children.end() && !found; ++tunIter) {
      NodeConn *conn = *tunIter;
      Node *tun = conn->m_n;
      if (!tun->IsTunnel(POSSTUNIN)) {
	LOG_FAIL("replacement for throw call");
	throw;
      }
      NodeConnVecIter iter = tun->m_children.begin();
      for(; iter != tun->m_children.end() && !found; ++iter) {
	NodeConn *childConn = *iter;
	if (childConn->m_num == subMat) {
	  if (!(childConn->m_n->IsTunnel(POSSTUNOUT))) {
	    found = true;
	    for(unsigned int execNum = 0; execNum < numExecs; ++execNum) {
	      const Size m = (*ms)[execNum];
	      const Size n = (*ns)[execNum];
	      
	      unsigned int numIters;
	      if (buildCache && !foundOne) {
		numIters = NumIters(bs, m, n);
		numItersVec.push_back(numIters);
	      }
	      else if (!foundOne) {
		numIters = NumIters(bs, m, n);
		if (numIters != numItersVec[execNum])
		  throw;
	      }
	      else {
		numIters = numItersVec[execNum];
	      }
	      
	      if (numIters) {
		GetSizes(subMat, numIters,
			 bs, parFactor,
			 m, n,
			 m_msizes[subMat], m_nsizes[subMat]);
	      }
	    }
	    foundOne = true;

	  }
	}
      }
    }
  }
  if (!foundOne) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
}
#else
void SplitSingleIter::BuildSizes(bool buildCache, vector<int> &numItersVec, unsigned int parFactor)
{
  if (m_tunType != SETTUNIN)
    return;
  if (m_pset && !m_pset->IsReal())
    return;
  if (!m_sizes) {
    LOG_FAIL("replacement for throw call");
    throw;
  }

  const Size bs = GetMyLoop()->GetBS();
  Dim numDims = InputNumDims(0);

  
  unsigned int numExecs = InputLen(0,0)->NumSizes();

  if (buildCache) {
    numItersVec.reserve(numExecs);
  }
  else {
    if (numExecs != numItersVec.size())
      LOG_FAIL("replacement for throw call");
  }


  for (Dim dim = 0; dim < numDims; ++dim) {
    const Sizes *sizes = InputLen(0,dim);
    unsigned int length = sizes->NumSizes();
    if (length != numExecs) {
      LOG_FAIL("replacement for throw call");
    }
    for(unsigned int execNum = 0; execNum < numExecs; ++execNum) {
      unsigned int numIters;
      if (buildCache && dim == 0) {
	numIters = NumIters(execNum);
	numItersVec.push_back(numIters);
      }
      else if (!buildCache || dim != 0) {
	numIters = numItersVec[execNum];
      }
      else {
	numIters = NumIters(execNum);
	if (numIters != numItersVec[execNum])
	  LOG_FAIL("replacement for throw call"); 
      }
      if (numIters) {
	const Size len = (*sizes)[execNum];
	if (dim < m_partDim) {
	  m_sizes[dim].AddRepeatedSizes(len, numIters, parFactor);
	}
	else if (dim == m_partDim) {
	  m_sizes[dim].AddSizesWithLimit(0,bs,len, parFactor);
	  m_sizes[dim+1].AddMidSizes(bs, len, parFactor);
	  m_sizes[dim+2].AddSizesWithLimit(len,-1*bs,0, parFactor);
	}
	else {
	  m_sizes[dim+2].AddRepeatedSizes(len, numIters, parFactor);
	}
      }
    }
  }
}
#endif
