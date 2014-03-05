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

#include "split.h"
#include "combine.h"
#include "distributions.h"
#include <cmath>

Split::Split() : LoopTunnel(LASTTUNNEL), 
#if TWOD
m_dir(LASTPARTDIR)
#else
m_partDim(99)
#endif
{
  m_addDir = false;
}

#if TWOD
Split::Split(PartDir dir, PossTunType type, bool isControl) 
  : LoopTunnel(type), m_dir(dir), m_isControlTun(isControl)
{
  m_addDir = false;
}
#else
Split::Split(unsigned int partDim, PossTunType type, bool isControl) 
  : LoopTunnel(type), m_partDim(partDim), m_isControlTun(isControl)
{
  m_addDir = false;
}
#endif

Split::~Split()
{
#if TWOD
  if (m_msizes) {
    delete [] m_msizes;
    m_msizes = NULL;
    delete [] m_mlsizes;
    m_mlsizes = NULL;
    delete [] m_nsizes;
    m_nsizes = NULL;
    delete [] m_nlsizes;
    m_nlsizes = NULL;
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

Node* Split::BlankInst() 
{ 
#if TWOD
  return new Split(LASTPARTDIR,LASTTUNNEL,false);
#else
  return new Split(99,LASTTUNNEL,false);
#endif
}

Name Split::GetName(unsigned int num) const 
{
  return GetName(num, GetLoopType());
}

Name Split::GetName(unsigned int num, LoopType type) const 
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
	throw;
      }
  }
#else
  if (m_addDir)
    throw;
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
        else 
          throw;
        return name;
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
        else
          throw;
        return name;
      default:
        cout << "bad dir\n";
        throw;
        break;
      }
#else
    if (m_partDim >= InputNumDims(0))
      throw;
    if (num > 2) 
      throw;
    else
      name.m_name += num;      
#endif
    return name;
  }
  else
    return name;
}

void Split::Prop()
{
  if (!IsValidCost(m_cost)) {
    for (unsigned int i = 0; i < m_inputs.size(); ++i)
      Input(i)->Prop();
    m_cost = ZERO;
    if (GetMyLoop()->GetBS() == ZERO)
      throw;
    if (m_tunType == POSSTUNIN) {
      if (Input(0)->GetNodeClass() != Split::GetClass())
        throw;
      Node *child = NULL;
      unsigned int size = m_children.size();  
      for (unsigned int i = 0; i < size; ++i) {
#if TWOD
        if (m_children[i]->m_num == GetNumElems(m_dir)) {
#else
        if (m_children[i]->m_num == 3) {
#endif
          if (child)
            throw;
          else
            child = m_children[i]->m_n;
        }
      }
      if (!child) {
        cout << GetNameStr(0) << endl;
        //m_poss->ForcePrint();
        throw;
      }
      else if (child->GetNodeClass() != Combine::GetClass())
        throw;
      else {
        Combine* com = (Combine*)child;
#if TWOD
        if (com->m_dir != m_dir)
#else
	if (com->m_partDim != m_partDim)
#endif
          throw;
      }
    }
  }
}

#if TWOD
const Sizes* Split::GetM(unsigned int num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0)
        throw;
      return GetInputM(0);
    case (POSSTUNIN):
      if (num < GetNumElems(m_dir)) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_msizes)
          throw;
        return &(input->m_msizes[num]);
      }
      else if (num == GetNumElems(m_dir)) {
        return GetInputM(0);
      }
      else
        throw;
    default:
      throw;
    }
}

const Sizes* Split::GetN(unsigned int num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0)
        throw;
      return GetInputN(0);
    case (POSSTUNIN):
      if (num < GetNumElems(m_dir)) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_nsizes)
          throw;
        return &(input->m_nsizes[num]);
      }
      else if (num == GetNumElems(m_dir)) {
        return GetInputN(0);
      }
      else
        throw;
    default:
      throw;
    }
}

const Sizes* Split::LocalM(unsigned int num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0)
        throw;
      return InputLocalM(0);
    case (POSSTUNIN):
      if (num < GetNumElems(m_dir)) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_mlsizes)
          throw;
        return input->m_mlsizes+num;
      }
      else if (num == GetNumElems(m_dir)) {
        return InputLocalM(0);
      }
      else
        throw;
    default:
      throw;
    }
}

const Sizes* Split::LocalN(unsigned int num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0)
        throw;
      return InputLocalN(0);
    case (POSSTUNIN):
      if (num < GetNumElems(m_dir)) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_nlsizes)
          throw;
        return input->m_nlsizes+num;
      }
      else if (num == GetNumElems(m_dir)) {
        return InputLocalN(0);
      }
      else
        throw;
    default:
      throw;
    }
}


void Split::GetSizes(unsigned int num, unsigned int numIters,
		     Size bs, unsigned int parFactor,
		     Size m, Size n,
		     Sizes &ms, Sizes &ns)
{
  if (m_tunType != SETTUNIN)
    throw;
  //  if (GetLoopType() != ELEMLOOP) {
  //    cout << "check that the m and n values are actually local sizes\n";
  //    cout << "also use currlm = currm...\n";
  //    throw;
  //  }
  switch(m_dir)
    {
    case (PARTDOWN):
      if (num == 0) {
	ms.AddSizesWithLimit(0,bs,m, parFactor);
	ns.AddRepeatedSizes(n,numIters, parFactor);
      }
      else if (num ==1) {
	if (ceil(m/bs) != numIters)
	  throw;
	ms.AddMidSizes(bs, m, parFactor);
	ns.AddRepeatedSizes(n,numIters, parFactor);
      }
      else if (num == 2) {
	ms.AddSizesWithLimit(m,-1*bs,0, parFactor);
	ns.AddRepeatedSizes(n,numIters, parFactor);
      }
      else 
        throw;
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
      else 
        throw;
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
      else
        throw;
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
      else
        throw;
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
      else
        throw;
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
      else
        throw;
      break;
    default:
      throw;
      break;
    }
}
#else
 const unsigned int Split::NumDims(unsigned int num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0)
        throw;
      return InputNumDims(0);
    case (POSSTUNIN):
      if (num > 3)
	throw;
      return InputNumDims(0);
    default:
      throw;
    }
}

const Sizes* Split::Len(unsigned int num, unsigned int dim) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0)
        throw;
      return InputLen(0,dim);
    case (POSSTUNIN):
      if (num < 3) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_sizes)
          throw;
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
      else
        throw;
    default:
      throw;
    }
}

   const Sizes* Split::LocalLen(unsigned int num, unsigned int dim) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
    case (SETTUNOUT):
      if (num > 0)
        throw;
      return InputLocalLen(0,dim);
    case (POSSTUNIN):
      if (num < 3) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_lsizes)
          throw;
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
      else
        throw;
    default:
      throw;
    }
}
#endif

#if TWOD
PossTunnel* Split::GetSetTunnel()
{
  Split *tun;
  if (m_tunType == POSSTUNIN)
    tun = new Split(m_dir, SETTUNIN);
  else if (m_tunType == POSSTUNOUT)
    tun = new Split(m_dir, SETTUNOUT);
  else
    throw;
  tun->CopyTunnelInfo(this);
  return tun;
}
#else
PossTunnel* Split::GetSetTunnel()
{
  Split *tun;
  if (m_tunType == POSSTUNIN)
    tun = new Split(m_partDim, SETTUNIN);
  else if (m_tunType == POSSTUNOUT)
    tun = new Split(m_partDim, SETTUNOUT);
  else
    throw;
  tun->CopyTunnelInfo(this);
  return tun;
}
#endif

void Split::PrintCode(IndStream &out)
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
	if (!child->IsPossTunnel(POSSTUNOUT)) {
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
	    else
	      throw;
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
	    else
	      throw;
	    break;
	  default:
	    throw;
	  }
	  string loopLevel = out.LoopLevel();
	  *out << "( " << part << ", idx" << loopLevel << ", bs" << loopLevel
	       << ", &" << inName << ", &" << GetNameStr(num) << " );\n";
	}
      }
    }
  }
#else
  *out << "need split print code\n";
  throw;
#endif
}

void Split::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  LoopTunnel::Duplicate(orig, shallow, possMerging);
  const Split *split = (Split*)orig;
#if TWOD
  m_dir = split->m_dir;
#else
  m_partDim = split->m_partDim;
#endif
  m_isControlTun = split->m_isControlTun;
  m_addDir = split->m_addDir;
}

NodeType Split::GetType() const
{
#if TWOD
  return "Split " + PartDirToStr(m_dir)  + "( " + PossTunnel::GetType() + " )";
#else
  string tmp = "Split";
  tmp += m_partDim;
  return tmp  + "( " + PossTunnel::GetType() + " )";
#endif
}

void Split::SanityCheck()
{
  LoopTunnel::SanityCheck();
  
#if TWOD
  if ((m_dir == PARTDOWN || m_dir == PARTUPWARD)
      && (GetUpStat(TL) != GetUpStat(TR) || GetUpStat(BL) != GetUpStat(BR))) 
    {
      cout << "bad statuses\n";
      throw;
    }
  else if ((m_dir == PARTRIGHT ||m_dir == PARTLEFT)
           && (GetUpStat(TL) != GetUpStat(BL) || GetUpStat(TR) != GetUpStat(BR)))
    {
      cout << "bad statuses\n";
      throw;
    }
#else
  throw;
#endif  
  if (m_tunType == POSSTUNIN) {
    for(unsigned int i = 0; i < m_inputs.size(); ++i) {
      if (Input(i)->GetNodeClass() != Split::GetClass()) {
        throw;
      }
    }
  }
  else if (m_tunType == SETTUNIN) {
    if (!m_pset->IsLoop())
      throw;
    if (m_inputs.size() != 1) {
      cout << "split has wrong number of inputs\n";
      throw;
    }
    for (unsigned int i = 0; i < m_children.size(); ++i) {
      if (Child(i)->GetNodeClass() != Split::GetClass()) {
        throw;
      }
    }
  }
  else {
    cout << "bad tunType\n";
    throw;
  }
}


unsigned int Split::NumOutputs() const 
{
#if TWOD
  return GetNumElems(m_dir)+1;
#else
  return 4;
#endif
}

#if TWOD
bool Split::QuadInUse(Quad quad, bool atEnd) const
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
      bool check = false;
      unsigned int num = (*iter)->m_num;
      switch (m_dir) {
      case (PARTDOWN):
	if ((quad == TL || quad == TR) && (!atEnd && num == 0 || atEnd && (num == 0 || num == 1)))
	  check = true;
	else if ((quad == BL || quad == BR) && (!atEnd && (num == 1 || num == 2) || atEnd && (num == 2)))
	  check = true;
	else
	  check = false;
	break;
      case (PARTUPWARD):
	if ((quad == TL || quad == TR) && (!atEnd && (num == 0 || num == 1) || atEnd && (num == 0)))
	  check = true;
	else if ((quad == BL || quad == BR) && (!atEnd && num == 2 || atEnd && (num == 1 || num == 2)))
	  check = true;
	else
	  check = false;
	break;
      case (PARTRIGHT):
	if ((quad == TL || quad == BL) && (!atEnd && num == 0 || atEnd && (num == 0 || num == 1)))
	  check = true;
	else if ((quad == TR || quad == BR) && (!atEnd && (num == 1 || num == 2) || atEnd && num == 2))
	  check = true;
	else
	  check = false;
	break;
      case (PARTLEFT):
	if ((quad == TL || quad == BL) && (!atEnd && (num == 0 || num == 1) || atEnd && num == 0))
	  check = true;
	else if ((quad == TR || quad == BR) && (!atEnd && num == 2 || atEnd && (num == 1 || num == 2)))
	  check = true;
	else
	  check = false;
	break;
      case (PARTDIAG):
	if (quad == TL && (!atEnd && num == 0 || atEnd && (num == 0 || num == 1 || num == 3 || num == 4)))
	  check = true;
	else if (quad == BL && (!atEnd && (num == 1 || num == 2) || atEnd && (num == 2 || num == 5)))
	  check = true;
	else if (quad == TR && (!atEnd && (num == 3 || num == 6) || atEnd && (num == 6 || num == 7)))
	  check = true;
	else if (quad == BR && (!atEnd && (num == 4 || num == 5 || num == 7 || num == 8)) || atEnd && num == 8)
	  check = true;
	else
	  check = false;
	break;
      case (PARTDIAGBACK):
	if (quad == TL && (!atEnd && (num == 0 || num == 1 || num == 3 || num == 4) || atEnd && num == 0))
	  check = true;
	else if (quad == BL && (!atEnd && (num == 2 || num == 5) || atEnd && (num == 1 || num == 2)))
	  check = true;
	else if (quad == TR && (!atEnd && (num == 6 || num == 7) || atEnd && (num == 3 || num == 6)))
	  check = true;
	else if (quad == BR && (!atEnd && (num == 8) || atEnd && (num == 4 || num == 5 || num == 7 || num == 8)))
	  check = true;
	else
	  check = false;
	break;
      default:
	throw;
      }
      if (check && !(*iter)->m_n->IsPossTunnel(POSSTUNOUT))
        return true;
    }
    return false;
  }
  else
    throw;
}
#else
bool Split::QuadInUse(Quad quad, bool atEnd) const
{
  throw;
}
#endif

void Split::PrintVarDeclarations(IndStream &out) const
{
  if (m_tunType != SETTUNIN)
    return;
  //BAMTODO: eventually predefine each obj_t matrix to be used ahead of the loop
  return;
#if TWOD
  LoopType type = GetLoopType();
  if (type == ELEMLOOP) {
    
  }
  else if (type == BLISLOOP) {
    string var = GetName(0).str();
    out.Indent();
    if (m_dir == PARTDOWN || m_dir == PARTUPWARD) 
      *out << "obj_t " << var << "T," << var << "B, " 
      << var << "0, " << var << "1, " << var << "2\n";
    else if (m_dir == PARTRIGHT || m_dir == PARTLEFT)
      *out << "obj_t " << var << "L," << var << "R, " 
      << var << "0, " << var << "1, " << var << "2\n";
    else if (m_dir == PARTDIAG || m_dir == PARTDIAGBACK)
      *out << "obj_t " << var << "TL," << var << "TR, " 
      << var << "BL," << var << "BR, \n" 
      << out.Tabs(1)
      << var << "0, " << var << "1, " << var << "2, "
      << var << "3, " << var << "4, " << var << "5, "
      << var << "6, " << var << "7, " << var << "8\n";
  }
  else
    throw;
#else
  cout << "need print var code\n";
  throw;
#endif
}

Combine* Split::CreateMatchingCombine(int numArgs, ...)
{
  int numComIns = NumOutputs();
  int j = 0;
#if TWOD
  Combine *com = new Combine(m_dir, POSSTUNOUT);
#else
  Combine *com = new Combine(m_partDim, POSSTUNOUT);
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
    if (j >= numComIns)
      throw;
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
unsigned int Split::NumIters(Size bs, Size m, Size n) const
{
  if (!bs)
    throw;
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
      throw;
  }
  throw;
}
#else
unsigned int Split::NumIters(Size bs, Size size) const
{
  if (!bs)
    throw;
  return ceil(((double)size)/bs);
}
#endif

#if TWOD
unsigned int Split::NumIters(unsigned int iterNum) const
{
  Size bs = GetMyLoop()->GetBS();
  const bool BLISLoop = (GetMyLoop()->GetType() == BLISLOOP);
  const Sizes *ms = (BLISLoop ? InputLocalM(0) : GetInputM(0));
  const Sizes *ns = (BLISLoop ? InputLocalN(0) : GetInputN(0));
  if (!ms || !ns) {
    if (Input(0)->m_flags & BUILDFLAG)
      cout << "has built\n";
    else
      cout << "hasn't built\n";
    cout << Input(0)->GetNodeClass() << endl;
    cout << Input(0) << endl;
    throw;
  }
  const Size m = (*ms)[iterNum];
  const Size n = (*ns)[iterNum];
  return NumIters(bs, m, n);
}
#else
unsigned int Split::NumIters(unsigned int iterNum) const
{
  Size bs = GetMyLoop()->GetBS();
  const Sizes *sizes = InputLen(0,m_partDim);
  if (!sizes) {
    if (Input(0)->m_flags & BUILDFLAG)
      cout << "has built\n";
    else
      cout << "hasn't built\n";
    cout << Input(0)->GetNodeClass() << endl;
    cout << Input(0) << endl;
    throw;
  }
  const Size m = (*sizes)[iterNum];
  return NumIters(bs,m);
}
#endif

void Split::FlattenCore(ofstream &out) const
{
  LoopTunnel::FlattenCore(out);
#if TWOD
  WRITE(m_dir);
#else
  WRITE(m_partDim);
#endif
  WRITE(m_isControlTun);
  WRITE(m_addDir);
}


void Split::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  LoopTunnel::UnflattenCore(in,info);
#if TWOD
  READ(m_dir);
#else
  READ(m_partDim);
#endif
  READ(m_isControlTun);
  READ(m_addDir);
}

unsigned int Split::NumberOfLoopExecs() const
{
  if (!m_isControlTun)
    throw;
#if TWOD
  //  return min(GetInputM(0)->NumSizes(), GetInputN(0)->NumSizes());
  unsigned int one = GetInputM(0)->NumSizes();
  unsigned int two = GetInputN(0)->NumSizes();
  if (one != two) {
    cout << one << " vs. " << two << endl;
    cout.flush();
    throw;
  }
  return min(one,two);
#else
  return InputLen(0,0)->NumSizes();
#endif
}

#if TWOD
void Split::StartFillingSizes()
{
  if (m_msizes)
    throw;
  if (m_tunType != SETTUNIN)
    return;
  unsigned int numElems = GetNumElems(m_dir);
  m_msizes = new Sizes[numElems];
  m_nsizes = new Sizes[numElems];
  m_mlsizes = new Sizes[numElems];
  m_nlsizes = new Sizes[numElems];
}
#else
void Split::StartFillingSizes()
{
  if (m_sizes)
    throw;
  if (m_tunType != SETTUNIN)
    return;
  unsigned int numDims = InputNumDims(0);
  //Num dims of sizes, but the m_partDim dimension has
  // 3 outputs
  m_sizes = new Sizes[numDims+2];
  m_lsizes = new Sizes[numDims+2];
}
#endif

void Split::ClearSizeCache()
{
#if TWOD
  if (!m_msizes)
    return;
  delete [] m_msizes;
  m_msizes = NULL;
  delete [] m_mlsizes;
  m_mlsizes = NULL;
  delete [] m_nsizes;
  m_nsizes = NULL;
  delete [] m_nlsizes;
  m_nlsizes = NULL;
#else
  if (!m_sizes)
    return;
  delete  [] m_sizes;
  m_sizes = NULL;
  delete [] m_lsizes;
  m_lsizes = NULL;
#endif
}

#if TWOD
void Split::AppendSizes(unsigned int execNum, unsigned int numIters, unsigned int parFactor)
{
  if (m_tunType != SETTUNIN)
    return;
  if (!m_msizes)
    throw;
  const LoopType loopType = GetMyLoop()->GetType();
  const bool BLISLoop = loopType == BLISLOOP;
  const Sizes *ms = (BLISLoop ? InputLocalM(0) : GetInputM(0));
  const Sizes *ns = (BLISLoop ? InputLocalN(0) : GetInputN(0));
  
  if (!ms || !ns) {
    if (Input(0)->m_flags & BUILDFLAG)
      cout << "has built\n";
    else
      cout << "hasn't built\n";

    cout << Input(0)->GetNodeClass() << endl;
    cout << Input(0) << endl;
    throw;
  }

  unsigned int length = ms->NumSizes();
  unsigned int length2 = ns->NumSizes();


  if (length <= execNum || length2 <= execNum)
    {
      throw;
    }
  const Size m = (*ms)[execNum];
  const Size n = (*ns)[execNum];
  const Size bs = GetMyLoop()->GetBS();
  const unsigned int numElems = GetNumElems(m_dir);

  if (NumIters(bs, m, n) != numIters) {
    InputLocalN(0)->Print();
    cout << endl;
    GetInputN(0)->Print();
    cout << endl;
    cout << NumIters(bs, m, n) << " vs. " << numIters << endl;
    for (unsigned int i = m_children.size(); i < m_children.size(); ++i) {
      cout << m_children[i]->m_n->GetNodeClass() << endl;
    }
    throw;
  }
  for (unsigned int subMat = 0; subMat < numElems; ++subMat) {
    GetSizes(subMat, numIters,
             bs, parFactor,
             m, n,
             m_msizes[subMat], m_nsizes[subMat]);
    /*
    if (m_msizes[subMat].NumSizes() != m_nsizes[subMat].NumSizes()) {
    cout << m_msizes[subMat].NumSizes() << endl;
    cout << m_nsizes[subMat].NumSizes() << endl;
      throw;
      }*/
  }
}
#else
void Split::AppendSizes(unsigned int execNum, unsigned int numIters, unsigned int parFactor)
{
  if (m_tunType != SETTUNIN)
    return;
  if (!m_sizes)
    throw;
  const Size bs = GetMyLoop()->GetBS();
  unsigned int numDims = InputNumDims(0);
  for (unsigned int dim; dim < numDims; ++dim) {
    const Sizes *sizes = InputLen(0,dim);
    unsigned int length = sizes->NumSizes();
    if (length <= execNum) {
      throw;
    }
    const Size len = (*sizes)[execNum];

    if (dim < m_partDim) {
      m_sizes[dim].AddRepeatedSizes(len, numIters, parFactor);
    }
    else if (dim == m_partDim) {
      if (NumIters(bs, len) != numIters)  {
	sizes->Print();
	cout << endl;
	cout << NumIters(bs, len) << " vs. " << numIters << endl;
	
      }
      m_sizes[dim].AddSizesWithLimit(0,bs,len, parFactor);
      m_sizes[dim+1].AddMidSizes(bs, len, parFactor);
      m_sizes[dim+2].AddSizesWithLimit(len,-1*bs,0, parFactor);
    }
    else {
      m_sizes[dim+2].AddRepeatedSizes(len, numIters, parFactor);
    }
  }
}
#endif

#if TWOD
void Split::UpdateLocalSizes()
{
  const unsigned int numElems = GetNumElems(m_dir);

  if (loopType == ELEMLOOP) {
#if DOELEM
    const DistType t = InputDistType(0);

    for (unsigned int subMat = 0; subMat < numElems; ++subMat) {
      GetLocalSizes(t, m_msizes+subMat, m_nsizes+subMat, m_mlsizes[subMat], m_nlsizes[subMat]);
    }
#endif
  }
  else {
    for (unsigned int subMat = 0; subMat < numElems; ++subMat) {
      m_mlsizes[subMat] = m_msizes[subMat];
      m_nlsizes[subMat] = m_nsizes[subMat];
    }    
  }
}
#else
void Split::UpdateLocalSizes()
{
  unsigned int numDims = InputNumDims(0);
  const DistType t = InputDistType(0);
  for (unsigned int dim = 0; dim < numDims; ++ dim) {
    GetLocalSizes(t, m_sizes+dim, m_lsizes+dim);
  }
}
#endif
