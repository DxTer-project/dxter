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

#include "splitUnrolled.h"
#include "combineUnrolled.h"
#include "elemRedist.h"
#include <cmath>
#include "LLDLA.h"

SplitUnrolled::SplitUnrolled() 
  : SplitBase(),
    m_unrollFactor(0)
{
}

#if TWOD
SplitUnrolled::SplitUnrolled(PartDir dir, unsigned int unrollFactor, TunType type, bool isControl) 
  : SplitBase(dir, type, isControl), m_unrollFactor(unrollFactor)
{
}
#else
SplitUnrolled::SplitUnrolled(unsigned int partDim, unsigned int unrollFactor, TunType type, bool isControl) 
  : SplitBase(partDim, type, isControl), m_unrollFactor(unrollFactor)
{
}
#endif


Node* SplitUnrolled::BlankInst() 
{ 
#if TWOD
  return new SplitUnrolled(LASTPARTDIR,0,LASTTUNNEL,false);
#else
  return new SplitUnrolled(99,0,LASTTUNNEL,false);
#endif
}

Name SplitUnrolled::GetName(ConnNum num) const 
{
  return GetName(num, GetLoopType());
}

Name SplitUnrolled::GetName(ConnNum num, LoopType type) const 
{
  Name name;
  if (type == BLISLOOP) {
    if (m_tunType == SETTUNOUT)
      if (m_pset->IsReal())
	return ((LoopTunnel*)Input(0))->GetOrigName();
      else
	return ((LoopTunnel*)(GetRealTunnel()->Input(0)))->GetOrigName();
    else if (m_tunType == POSSTUNOUT)
      return ((LoopTunnel*)Input(m_inputs.size()-1))->GetOrigName();
    else if (m_tunType == SETTUNIN) {
      name.m_name = GetInputName(0).str();
      if (name.m_name[name.m_name.length()-1] != '_')
        name.m_name = name.m_name + "_";
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
    if (num < m_unrollFactor)
      name.m_name += std::to_string((long long int) num);
    else if (num != m_unrollFactor)
      throw;
    return name;
  }
  else
    return name;
}

void SplitUnrolled::Prop()
{
  if (!IsValidCost(m_cost)) {
    Tunnel::Prop();
    if (!m_unrollFactor)
      throw;
    m_cost = ZERO;
  }
}

#if TWOD
const Sizes* SplitUnrolled::GetM(ConnNum num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
      if (num > 0)
        throw;
      return GetInputM(0);
    case (SETTUNOUT):
      return GetRealTunnel()->GetInputM(0);
    case (POSSTUNIN):
      if (num < m_unrollFactor) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_msizes)
          throw;
        return input->m_msizes;
      }
      else if (num == m_unrollFactor) {
        return GetInputM(0);
      }
      else
        throw;
    default:
      throw;
    }
}

const Sizes* SplitUnrolled::GetN(ConnNum num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
      if (num > 0)
        throw;
      return GetInputN(0);
    case (SETTUNOUT):
      return GetRealTunnel()->GetInputN(0);
    case (POSSTUNIN):
      if (num < m_unrollFactor) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_nsizes)
          throw;
        return input->m_nsizes+num;
      }
      else if (num == m_unrollFactor) {
        return GetInputN(0);
      }
      else
        throw;
    default:
      throw;
    }
}

#if DODM
const Sizes* SplitUnrolled::LocalM(ConnNum num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
      if (num > 0)
        throw;
      return InputLocalM(0);
    case (SETTUNOUT):
      GetRealTunnel()->InputLocalM(0);
    case (POSSTUNIN):
      if (num < m_unrollFactor) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_mlsizes)
          throw;
        return input->m_mlsizes;
      }
      else if (num == m_unrollFactor) {
        return InputLocalM(0);
      }
      else
        throw;
    default:
      throw;
    }
}

const Sizes* SplitUnrolled::LocalN(ConnNum num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
      if (num > 0)
        throw;
      return InputLocalN(0);
    case (SETTUNOUT):
      return GetRealTunnel()->InputLocalN(0);
    case (POSSTUNIN):
      if (num < m_unrollFactor) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_nlsizes)
          throw;
        return input->m_nlsizes;
      }
      else if (num == m_unrollFactor) {
        return InputLocalN(0);
      }
      else
        throw;
    default:
      throw;
    }
}
#endif



#else
 const Dim SplitUnrolled::NumDims(ConnNum num) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
      if (num > 0)
        throw;
      return InputNumDims(0);
    case (SETTUNOUT):
      GetRealTunnel()->InputNumDims(0);
    case (POSSTUNIN):
      if (num > m_unrollFactor)
	throw;
      return InputNumDims(0);
    default:
      throw;
    }
}

const Sizes* SplitUnrolled::Len(ConnNum num, Dim dim) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
      if (num > 0)
        throw;
      return InputLen(0,dim);
    case (SETTUNOUT):
      return GetRealTunnel()->InputLen(0,dim);
    case (POSSTUNIN):
      if (num < m_unrollFactor) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_sizes)
          throw;
	return &(input->m_sizes[dim]);
      }
      else if (num == m_unrollFactor) {
        return InputLen(0,dim);
      }
      
      else
        throw;
    default:
      throw;
    }
}

const Sizes* SplitUnrolled::LocalLen(ConnNum num, Dim dim) const
{
  switch(m_tunType) 
    {
    case (SETTUNIN):
    case (POSSTUNOUT):
      if (num > 0)
        throw;
      return InputLocalLen(0,dim);
    case (SETTUNOUT):
      GetRealTunnel()->InputLocalLen(0,dim);
    case (POSSTUNIN):
      if (num < m_unrollFactor) {
        const LoopTunnel *input = (LoopTunnel*)Input(0);
        if (!input->m_lsizes)
          throw;
	return &(input->m_lsizes[dim]);
      }
      else if (num == m_unrollFactor) {
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
Tunnel* SplitUnrolled::GetSetTunnel()
{
  SplitUnrolled *tun;
  if (m_tunType == POSSTUNIN)
    tun = new SplitUnrolled(m_dir, m_unrollFactor, SETTUNIN);
  else if (m_tunType == POSSTUNOUT)
    tun = new SplitUnrolled(m_dir, m_unrollFactor, SETTUNOUT);
  else
    throw;
  tun->CopyTunnelInfo(this);
  return tun;
}
#else
Tunnel* SplitUnrolled::GetSetTunnel()
{
  SplitUnrolled *tun;
  if (m_tunType == POSSTUNIN)
    tun = new SplitUnrolled(m_partDim, m_unrollFactor, SETTUNIN);
  else if (m_tunType == POSSTUNOUT)
    tun = new SplitUnrolled(m_partDim, m_unrollFactor, SETTUNOUT);
  else
    throw;
  tun->CopyTunnelInfo(this);
  return tun;
}
#endif

void SplitUnrolled::PrintCode(IndStream &out)
{
  if (m_tunType != POSSTUNIN)
    return;
#if TWOD
  LoopType loopType = GetLoopType();
  if (loopType == ELEMLOOP) {
    throw;
  }
  else if (loopType == BLISLOOP) {
    throw;
  }
#else
  *out << "need split print code\n";
  throw;
#endif
}

void SplitUnrolled::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  LoopTunnel::Duplicate(orig, shallow, possMerging);
  const SplitUnrolled *split = (SplitUnrolled*)orig;
#if TWOD
  m_dir = split->m_dir;
#else
  m_partDim = split->m_partDim;
#endif
  m_unrollFactor = split->m_unrollFactor;
  m_isControlTun = split->m_isControlTun;
}

NodeType SplitUnrolled::GetType() const
{
#if TWOD
  return "SplitUnrolled" + PartDirToStr(m_dir) + std::to_string((long long int) m_unrollFactor)
    + "( " + Tunnel::GetType() + " )";
#else
  string tmp = "SplitUnrolled";
  tmp += m_partDim;
  tmp += std::to_string((long long int) m_unrollFactor);
  return tmp  + "( " + Tunnel::GetType() + " )";
#endif
}

#if TWOD
bool SplitUnrolled::QuadInUse(Quad quad, bool atEnd) const
{
  throw;
}
#else
bool SplitUnrolled::QuadInUse(Quad quad, bool atEnd) const
{
  throw;
}
#endif

void SplitUnrolled::PrintVarDeclarations(IndStream &out) const
{
#if DOLLDLA
  if (m_tunType != POSSTUNIN)
    throw;
  const string name = GetInputNameStr(0);
  throw;
  for (unsigned int i = 0; i < m_unrollFactor; ++i) {
    out.Indent();
    *out << LLDLAPartVarName(name, i) << " = " 
	 << name << ";\n";
  }
#endif
}

void SplitUnrolled::AddVariables(VarSet &set) const
{
#if DOLLDLA
  if (m_tunType != POSSTUNIN)
    return;
  if (GetLoopType() == LLDLALOOP) {
    const string name = GetInputNameStr(0);
    for(unsigned int i = 0; i < m_unrollFactor; ++i) {
      Var var(name, i);
      set.insert(var);
    }
  }
  else
    throw;
#endif
}

CombineUnrolled* SplitUnrolled::CreateMatchingCombine(int numArgs, ...)
{
  throw;
  int numComIns = NumOutputs();
  int j = 0;
#if TWOD
  CombineUnrolled *com = new CombineUnrolled(m_dir, m_unrollFactor, POSSTUNOUT);
#else
  CombineUnrolled *com = new CombineUnrolled(m_partDim, m_unrollFactor, POSSTUNOUT);
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
unsigned int SplitUnrolled::NumIters(Size bs, Size m, Size n) const
{
  if (!bs)
    throw;
  switch (m_dir) {
    case(PARTDOWN):
    case(PARTDIAG):
    case(PARTUPWARD):
    case(PARTDIAGBACK):
      {
	double tmp = ((double)m)/bs;
	if (fmod(tmp, m_unrollFactor))
	  throw;
	else
	  return tmp / m_unrollFactor;
      }
    case(PARTLEFT):
    case(PARTRIGHT):
      {
	double tmp = ((double)n)/bs;
	if (fmod(tmp, m_unrollFactor))
	  throw;
	else
	  return tmp / m_unrollFactor;
      }
    case(LASTPARTDIR):
      throw;
  }
  throw;
}
#else
unsigned int SplitUnrolled::NumIters(Size bs, Size size) const
{
  if (!bs)
    throw;
  double tmp = ((double)size)/bs;
  if (fmod(tmp, m_unrollFactor))
    throw;
  else
    return tmp / m_unrollFactor;
}
#endif

#if TWOD
unsigned int SplitUnrolled::NumIters(unsigned int iterNum) const
{
  Size bs = GetMyLoop()->GetBS();

  const Sizes *ms = GetInputM(0);
  const Sizes *ns = GetInputN(0);
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
unsigned int SplitUnrolled::NumIters(unsigned int iterNum) const
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

void SplitUnrolled::FlattenCore(ofstream &out) const
{
  SplitBase::FlattenCore(out);
  WRITE(m_unrollFactor);
}


void SplitUnrolled::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  SplitBase::UnflattenCore(in,info);
  READ(m_unrollFactor);
}

unsigned int SplitUnrolled::NumberOfLoopExecs() const
{
  if (!m_isControlTun)
    throw;
#if TWOD
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
void SplitUnrolled::StartFillingSizes()
{
  if (m_msizes)
    throw;
  if (m_tunType != SETTUNIN)
    return;
  if (!m_pset->IsReal())
    return;
  m_msizes = new Sizes;
  m_nsizes = new Sizes;
#if DODM
  m_mlsizes = new Sizes;
  m_nlsizes = new Sizes;
#endif
}
#else
void SplitUnrolled::StartFillingSizes()
{
  if (m_sizes)
    throw;
  if (m_tunType != SETTUNIN)
    return;
  if (!m_pset->IsReal())
    return;
  unsigned int numDims = InputNumDims(0);
  //Num dims of sizes, but the m_partDim dimension has
  // m_unrollFactor outputs
  m_sizes = new Sizes[numDims];
  m_lsizes = new Sizes[numDims];
}
#endif

void SplitUnrolled::ClearDataTypeCache()
{
#if TWOD
  if (!m_msizes)
    return;
  delete m_msizes;
  m_msizes = NULL;
  delete m_nsizes;
  m_nsizes = NULL;
#if DODM
  delete m_mlsizes;
  m_mlsizes = NULL;
  delete m_nlsizes;
  m_nlsizes = NULL;
#endif
#else
  if (!m_sizes)
    return;
  delete [] m_sizes;
  m_sizes = NULL;
  delete [] m_lsizes;
  m_lsizes = NULL;
#endif
}

#if TWOD
void SplitUnrolled::AppendSizes(unsigned int execNum, unsigned int numIters, unsigned int parFactor)
{
  if (m_tunType != SETTUNIN)
    return;
  if (!m_pset->IsReal())
    return;
  if (!m_msizes)
    throw;
  const Sizes *ms = GetInputM(0);
  const Sizes *ns = GetInputN(0);
  
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

  if (NumIters(bs, m, n) != numIters) {
    GetInputN(0)->Print();
    cout << endl;
    cout << NumIters(bs, m, n) << " vs. " << numIters << endl;
    for (unsigned int i = m_children.size(); i < m_children.size(); ++i) {
      cout << m_children[i]->m_n->GetNodeClass() << endl;
    }
    throw;
  }
  m_msizes->AddRepeatedSizes(bs, numIters, 1);
}
#else
void SplitUnrolled::AppendSizes(unsigned int execNum, unsigned int numIters, unsigned int parFactor)
{
  if (m_tunType != SETTUNIN)
    return;
  if (!m_pset->IsReal())
    return;
  if (!m_sizes)
    throw;
  const Size bs = GetMyLoop()->GetBS();
  Dim numDims = InputNumDims(0);
  for (Dim dim = 0; dim < numDims; ++dim) {
    const Sizes *sizes = InputLen(0,dim);
    unsigned int length = sizes->NumSizes();
    if (length <= execNum) {
      throw;
    }
    const Size len = (*sizes)[execNum];
    
    if (dim == m_partDim) {
      if (NumIters(bs, len) != numIters)  {
	sizes->Print();
	cout << endl;
	cout << NumIters(bs, len) << " vs. " << numIters << endl;
	
      }
      m_sizes[dim].AddRepeatedSizes(bs, numIters, 1);
    }
    else {
      m_sizes[dim].AddRepeatedSizes(len, numIters, parFactor);
    }
  }
}
#endif

#if TWOD&&DODM
void SplitUnrolled::UpdateLocalSizes()
{
  throw;
}
#elif DOTENSORS
void SplitUnrolled::UpdateLocalSizes()
{
  throw;
}
#endif

#if DOLLDLA
string SplitUnrolled::LoopBound()
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
      throw;
      break;
    case (PARTUPWARD):
      throw;
      break;
    case (PARTLEFT):     
      throw;
      break;

    case (PARTDIAGBACK):
      throw;
      break;
    default:
      throw;
    }
}
#endif


void SplitUnrolled::PrintIncrementAtEndOfLoop(BSSize bs, IndStream &out) const
{
  if (m_tunType != POSSTUNIN)
    throw;
#if DOLLDLA
  throw;
  //need bs
  out.Indent(1);
  const DataTypeInfo &type = InputDataType(0);
  for (unsigned int i = 0; i < m_unrollFactor; ++i) {
    *out << LLDLAPartVarName(GetInputNameStr(0),i) << " += ";
    if (m_dir == PARTDOWN) {
      if (!IsUnitStride(type.m_rowStride))
	*out << type.m_rowStrideVar << " * ";
      *out << MU_VAR_NAME << ";\n";
    }
    else if (m_dir == PARTRIGHT) {
      if (!IsUnitStride(type.m_colStride))
	*out << type.m_colStrideVar << " * ";
      *out << MU_VAR_NAME << ";\n";
    }
    else
      throw;
  }
#endif
}
