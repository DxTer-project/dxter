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



#include "tempVarNode.h"
#include "elemRedist.h"
#include <cmath>
#include "loopSupport.h"

#define THROWWHENCANTMOVETEMPVARNODE 0

TempVarNode::TempVarNode()
#if DOELEM
: m_info(D_LASTDIST), m_mlsize(NULL), m_nlsize(NULL)
#elif DOTENSORS
: m_info(), m_lsizes(NULL), m_sumLens(NULL)
#endif
{
  
}

TempVarNode::TempVarNode(string name)
#if DOELEM
: m_info(D_LASTDIST), m_name(name)
#elif DOTENSORS
: m_info(), m_name(name)
#else
: m_name(name)
#endif
#if DOELEM
, m_mlsize(NULL), m_nlsize(NULL)
#elif DOTENSORS
, m_lsizes(NULL), m_sumLens(NULL)
#endif
{}

#if DODM
TempVarNode::TempVarNode(DistType dist)
: m_info(dist),
#if TWOD
m_mlsize(NULL), m_nlsize(NULL)
#elif DOTENSORS
m_lsizes(NULL), m_sumLens(NULL)
#endif
{}
#endif


#if DODM
TempVarNode::TempVarNode(DistType dist, string name)
:  m_info(dist), m_name(name),
#if TWOD
m_mlsize(NULL), m_nlsize(NULL)
#elif DOTENSORS
m_lsizes(NULL),
m_sumLens(NULL)
#endif
{}
#endif




#if DOTENSORS
TempVarNode::TempVarNode(DistType dist, EntryList sumDims)
: m_lsizes(NULL),
m_sumLens(NULL),
m_sumDims(sumDims)
{
  //update below, too
  Dim numSumDims = sumDims.size();
  if (!numSumDims) {
    cout << "!numSumDims\n";
    throw;
  }
  DistType temp;
  temp.PrepForNumDims(dist.m_numDims+numSumDims);
  for (Dim dim = 0; dim < dist.m_numDims; ++dim)
    temp.m_dists[dim] = dist.m_dists[dim];
  EntryListIter iter = m_sumDims.begin();
  for (Dim dim = 0; iter != m_sumDims.end(); ++iter, ++dim) {
    temp.m_dists[dist.m_numDims + dim] = *iter;
  }
  temp.m_notReped = dist.m_notReped;
  m_info.SetDistAndClearPerm(temp);
}

TempVarNode::TempVarNode(DistType dist, EntryList  sumDims, string name)
:  m_name(name),
m_lsizes(NULL),
m_sumLens(NULL),
m_sumDims(sumDims)
{
  //update above, too
  Dim numSumDims = sumDims.size();
  if (!numSumDims) {
    cout << "!numSumDims 2\n";
    throw;
  }
  DistType temp;
  temp.PrepForNumDims(dist.m_numDims+numSumDims);
  for (Dim dim = 0; dim < dist.m_numDims; ++dim)
    temp.m_dists[dim] = dist.m_dists[dim];
  EntryListIter iter = m_sumDims.begin();
  for (Dim dim = 0; iter != m_sumDims.end(); ++iter, ++dim) {
    temp.m_dists[dist.m_numDims + dim] = *iter;
  }
  temp.m_notReped = dist.m_notReped;
  m_info.SetDistAndClearPerm(temp);
}
#endif


NodeType TempVarNode::GetType() const
{
  if (m_inputs.size() != 1) {
    cout << "m_inputs.size() != 1\n";
    cout.flush();
    throw;
  }
  
  return GetName(0).m_name;
}

void TempVarNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig,shallow, possMerging);
  const TempVarNode *temp = (TempVarNode*)orig;
  m_info = temp->m_info;
  m_name = temp->m_name;
#if DOTENSORS
  m_sumDims = temp->m_sumDims;
  m_align = temp->m_align;
  m_alignModes = temp->m_alignModes;
  m_alignModesSrc = temp->m_alignModesSrc;  
#endif
}

TempVarNode::~TempVarNode()
{
#if TWOD
#if DODM
  if (m_mlsize) {
    delete m_mlsize;
    m_mlsize = NULL;
    delete m_nlsize;
    m_nlsize = NULL;
  }
#endif
#elif DOTENSORS
  if (m_lsizes) {
    delete [] m_lsizes;
    m_lsizes = NULL;
  }
#endif
}

void TempVarNode::PrintCode(IndStream &out)
{
  /*
   if (GetLayer() == SQ2LAYER || GetLayer() == SQ1LAYER) {
   string name = GetInputNameStr(0);
   out.Indent();
   *out << "bli_obj_create( BLIS_DOUBLE, bli_obj_length("
   << name << "), bli_obj_width(" << name << "), 0, 0, &"
   << GetNameStr(0) << " );\n";
   out.Indent();
   *out << "bli_copym(&" << name << ", &" << GetNameStr(0) << ");\n";
   }*/
#if DOTENSORS
  if (!m_align.empty()) {
    out.Indent();
    *out << GetNameStr(0) << ".AlignModesWith( "
	 << ModeArrayVarName(m_alignModes) << ", "
	 << m_align << ", "
	 << ModeArrayVarName(m_alignModesSrc) << " );\n";
  }
  out.Indent();
  *out << "tempShape = " << GetInputNameStr(0) << ".Shape();\n";
  EntryListIter iter = m_sumDims.begin();
  for (; iter != m_sumDims.end(); ++iter) {
    out.Indent();
    DistEntry entry = *iter;
    DimVec vec = entry.DistEntryDims();
    if (vec.empty())
      throw;
    if (vec.size() == 1) {
      *out << "tempShape.push_back( g.Shape()[" << vec[0] << "] );\n";
    }
    else {
      *out << "tempShape.push_back( ";
      for( int i = 0; i < vec.size(); ++i) {
        if (i)
          *out << " * ";
        *out << "g.Shape()[" << vec[i] << "]";
      }
      *out << " );\n";
    }
  }
  out.Indent();
  *out << GetNameStr(0) << ".ResizeTo( tempShape );\n";
#endif
}


void TempVarNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (m_inputs.size() != 1) {
      cout << "m_inputs.size() != 1\n";
      cout.flush();
      throw;
    }
    Input(0)->Prop();
#if DOTENSORS
    m_cost = -1;
    //MaxNumberOfLocalElements(0);
#else
    m_cost = ZERO;
#endif
  }
}

const DataTypeInfo& TempVarNode::DataType(ConnNum num) const
{
  return m_info;
}

#if TWOD
const Sizes* TempVarNode::GetM(ConnNum num) const
{
  if (num > 0)
    throw;
  return GetInputM(0);
}

const Sizes* TempVarNode::GetN(ConnNum num) const
{
  if (num > 0)
    throw;
  return GetInputN(0);
}

#if DODM
const Sizes* TempVarNode::LocalM(ConnNum num) const
{
  if (num > 0)
    throw;
  return m_mlsize;
}

const Sizes* TempVarNode::LocalN(ConnNum num) const
{
  if (num > 0)
    throw;
  return m_nlsize;
}
#endif //DODM

#elif DOTENSORS
const Dim TempVarNode::NumDims(ConnNum num) const
{
  return InputNumDims(0) + m_sumDims.size();
}

const Sizes* TempVarNode::Len(ConnNum num, Dim reqDim) const
{
  if (num > 0)
    throw;
  Dim dim = m_info.GetPerm().MapFinishToStart(reqDim);
  if (!m_sumDims.empty()) {
    Dim dims = m_info.GetDist().m_numDims-m_sumDims.size();
    if (dim >= dims)
      return m_sumLens + (dim - dims);
    else
      return InputLen(0, dim);
  }
  else
    return InputLen(0, dim);
}

const Sizes* TempVarNode::LocalLen(ConnNum num, Dim reqDim) const
{
  if (num > 0)
    throw;
  Dim dim = m_info.GetPerm().MapFinishToStart(reqDim);
  if (!m_sumDims.empty()) {
    Dim dims = m_info.GetDist().m_numDims-m_sumDims.size();
    if (dim >= dims)
      return &m_ones;
    else
      return &(m_lsizes[dim]);
  }
  else
    return &(m_lsizes[dim]);
}
#endif

Name TempVarNode::GetName(ConnNum num) const
{
  if (num > 0)
    throw;
  Name tmp;
  if (m_name.empty()) {
    tmp = GetInputName(0);
#if DOTENSORS
    
#else
    tmp.m_name = tmp.m_name + "temp";
#endif
  }
  else {
    tmp.m_name = m_name;
  }
#if DODM
  tmp.m_type = m_info.GetDist();
#if DOTENSORS
  tmp.m_permutation = m_info.GetPerm();
#endif
#endif
  return tmp;
}

void TempVarNode::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  throw; // m_info
  out << m_name << endl;
#if DOTENSORS
  throw; //m_sumDims
#endif
}

void TempVarNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLANode::UnflattenCore(in, info);
  throw; // m_info
  getline(in,m_name);
#if DOTENSORS
  throw; //m_sumDims
#endif
}

void TempVarNode::ClearDataTypeCache()
{
#if TWOD&&DODM
  if (!m_mlsize)
    return;
  delete m_mlsize;
  m_mlsize = NULL;
  delete m_nlsize;
  m_nlsize = NULL;
#elif DOTENSORS
  if (m_lsizes)
    return;
  delete [] m_lsizes;
  m_lsizes = NULL;
  m_ones.ClearSizes();
  delete [] m_sumLens;
  m_sumLens = NULL;
#endif
}

void TempVarNode::BuildDataTypeCache()
{
#if DOELEM
  if (m_mlsize)
    return;
  m_mlsize = new Sizes;
  m_nlsize = new Sizes;
  GetLocalSizes(m_info.m_dist, GetM(0), GetN(0), *m_mlsize, *m_nlsize);
#elif DOTENSORS
  
  if (m_lsizes)
    return;
  Dim numDims = m_info.GetDist().m_numDims-m_sumDims.size();
  m_lsizes = new Sizes[numDims];
  
  DistType type = m_info.GetEffectiveDist();
  
  for (Dim dim = 0; dim < numDims; ++dim)
    GetLocalSizes(type, dim, InputLen(0,dim), m_lsizes+dim);
  
  
  m_sumLens = new Sizes[m_sumDims.size()];
  m_ones.AddRepeatedSizes(1, InputLen(0,0)->NumSizes(), 1);
  EntryListIter iter = m_sumDims.begin();
  for(Dim dim = 0; iter != m_sumDims.end(); ++dim, ++iter) {
    DimVec vec = (*iter).DistEntryDims();
    unsigned int numProcs = 1;
    DimVecIter iter2 = vec.begin();
    for(; iter2 != vec.end(); ++iter2)
      numProcs *= GridLens[*iter2];
    m_sumLens[dim].AddRepeatedSizes(numProcs, InputLen(0,0)->NumSizes(), 1);
  }
  
#endif
}

#if DOTENSORS

bool HasDecendentForApplication(const Node *node)
{
  if (node->IsTunnel()) {
    const Tunnel *tun = (Tunnel*)node;
    if (tun->m_tunType == POSSTUNOUT || tun->m_tunType == SETTUNOUT)
      return false;
    if (tun->m_tunType == SETTUNIN) {
      if (tun->IsLoopTunnel()) {
        const LoopTunnel *out = ((LoopTunnel*)tun)->GetMatchingOutTun();
        if (out->m_children.empty()) {
          if (!tun->m_pset)
            throw;
          if (tun->m_pset->IsShadow()) {
#if THROWWHENCANTMOVETEMPVARNODE
            throw;
#else
	    return false;
#endif
	  }
          else if (!((RealPSet*)(tun->m_pset))->m_shadows.empty()) {
#if THROWWHENCANTMOVETEMPVARNODE
            cout << node->GetNameStr(0) << endl;
            throw;
#else
	    return false;
#endif
          }
          return true;
        }
        else
          return false;
      }
      else {
        if (!tun->m_pset)
          throw;
        if (tun->m_pset->IsShadow()) {
#if THROWWHENCANTMOVETEMPVARNODE
          if (HasDecendentForApplication(tun->GetRealTunnel()))  {
            cout << node->GetNameStr(0) << endl;
            throw;
	  }
#endif
	    return false;
        }
        //go through each poss tuns
        NodeConnVecConstIter iter = tun->m_children.begin();
        for(; iter != tun->m_children.end(); ++iter) {
          if (HasDecendentForApplication((*iter)->m_n)) {
            if (!((RealPSet*)(tun->m_pset))->m_shadows.empty()) {
#if THROWWHENCANTMOVETEMPVARNODE
            throw;
#else
	    return false;
#endif
	    }
            return true;
          }
        }
        return false;
      }
    }
    else if (tun->m_tunType == POSSTUNIN) {
      NodeConnVecConstIter iter = tun->m_children.begin();
      for( ; iter != tun->m_children.end(); ++iter) {
        const NodeConn *conn = *iter;
        if (conn->m_num == 0) {
          if (HasDecendentForApplication(conn->m_n))
            return true;
        }
      }
      return false;
    }
    else
      return false;
  }
  else
    return false;
}

bool MoveTempVarNodeIntoLoop::CanApply(const Node *node) const
{
  if (node->m_children.size() > 1) {
    return false;
  }
  return HasDecendentForApplication(node->Child(0));
}

void MoveIn(Node *node, Node *newSrc, ConnNum newSrcNum, TempVarNode *tempVarNode)
{
  if (node->IsTunnel()) {
    Tunnel *tun = (Tunnel*)node;
    if (tun->m_tunType == SETTUNIN) {
#if USESHADOWS
      if (tun->m_pset->IsShadow() || !((RealPSet*)(tun->m_pset))->m_shadows.empty())
        throw;
#endif
      if (tun->IsLoopTunnel()) {
        LoopTunnel *out = ((LoopTunnel*)tun)->GetMatchingOutTun();
        if (out->m_children.empty()) {
          if (!tun->m_pset)
            throw;
#if USESHADOWS
          if (tun->m_pset->IsShadow())
            throw;
          else if (!((RealPSet*)(tun->m_pset))->m_shadows.empty())
            throw;
          ((RealPSet*)(tun->m_pset))->DisconnectFromSetsForMergingRecord();
#endif
          LoopTunnel *newSrcTunIn = NULL;
          // If the tempvar goes into this loop, find if the src for the tempvar
          // also goes into the loop in the same way (i.e., split the same way or
          // not) and use that. Otherwise, create a similar input to the loop.
          NodeConnVecIter iter = newSrc->m_children.begin();
          for(; !newSrcTunIn && iter != newSrc->m_children.end(); ++iter) {
            if ((*iter)->m_num == newSrcNum) {
              if ((*iter)->m_n->IsTunnel(SETTUNIN)) {
                Tunnel *inTun = (Tunnel*)((*iter)->m_n);
                if (inTun->m_pset == tun->m_pset) {
                  if (!inTun->IsLoopTunnel())
                    throw;
                  if (tun->IsSplit()) {
                    if (inTun->IsSplit()) {
                      if (((SplitBase*)inTun)->m_partDim ==
                          ((SplitBase*)tun)->m_partDim) {
                        newSrcTunIn = (LoopTunnel*)inTun;
                      }
                    }
                  }
                  else if (!inTun->IsSplit()) {
                    newSrcTunIn = (LoopTunnel*)inTun;
                  }
                }
              }
            }
          }
          if (!newSrcTunIn) {
            //create new tunnel into the loop?
            throw;
          }
          if (tun->m_children.empty())
            throw;
          //the name in the tempVarNode doesn't have _part0_part1_ ... etc,
          // but if sets on one of the posses are shadows or have shadows, then
          // the input temp edges must have the same names, so use the previous name here
          string name = tun->Child(0)->GetName(tun->IsSplit() ? 1 : 0).m_name;
          for(unsigned int i = 0; i < tun->m_children.size(); ++i) {
            Node *tempVarPossTunIn = tun->Child(i);
            Node *newSrcPossTunIn = newSrcTunIn->Child(i);
            if (tun->IsSplit()) {
              TempVarNode *newTemp = new TempVarNode(tempVarNode->m_info.GetDist(), name);
              newTemp->AddInput(newSrcPossTunIn, 1);
              if (newTemp->GetNameStr(0) != tempVarPossTunIn->GetNameStr(1)) {
                cout << newTemp->GetNameStr(0) << endl;
                cout << tempVarPossTunIn->GetNameStr(1) << endl;
                throw;
              }
              tempVarPossTunIn->RedirectChildren(1, newTemp, 0);
	      if (newTemp->m_children.size() == 1) {
		((DLANode*)(newTemp->Child(0)))->AlignInfo(newTemp->m_align,
					     newTemp->m_alignModes,
					     newTemp->m_alignModesSrc);
	      }
              tempVarPossTunIn->m_poss->AddNode(newTemp);
              tempVarPossTunIn->m_poss->m_fullyExpanded = false;
            }
            else {
              TempVarNode  *newTemp = new TempVarNode(tempVarNode->m_info.GetDist(), name);
              newTemp->AddInput(newSrcPossTunIn, 0);
              if (newTemp->GetNameStr(0) != tempVarPossTunIn->GetNameStr(0)) {
                cout << newTemp->GetNameStr(1) << endl;
                cout << tempVarPossTunIn->GetNameStr(1) << endl;
                throw;
              }
              tempVarPossTunIn->RedirectChildren(0, newTemp, 0);
	      if (newTemp->m_children.size() == 1) {
		((DLANode*)(newTemp->Child(0)))->AlignInfo(newTemp->m_align,
					     newTemp->m_alignModes,
					     newTemp->m_alignModesSrc);
	      }
              tempVarPossTunIn->m_poss->AddNode(newTemp);
              tempVarPossTunIn->m_poss->m_fullyExpanded = false;
            }
          }
          LoopTunnel *setTunIn = (LoopTunnel*)tun;
          LoopTunnel *setTunOut = setTunIn->GetMatchingOutTun();
          if (setTunIn->IsSplit() && ((SplitBase*)setTunIn)->m_isControlTun) {
            bool found = false;
            TunVecIter inIter = tun->m_pset->m_inTuns.begin();
            for(; !found && inIter != tun->m_pset->m_inTuns.end(); ++inIter) {
              if (*inIter != tun && ((LoopTunnel*)(*inIter))->IsSplit()) {
                SplitBase *split = ((SplitBase*)(*inIter));
                split->m_isControlTun = true;
                for (int i = 0; i < split->m_children.size(); ++i)
                  ((SplitBase*)(split->Child(i)))->m_isControlTun = true;
                found = true;
              }
            }
          }
          
          
          while (!setTunIn->m_children.empty()) {
            Node *tempVarPossTunIn = setTunIn->Child(0);
            tempVarPossTunIn->RemoveAllChildren2Way();
            delete tempVarPossTunIn->m_inputs[0];
            tempVarPossTunIn->m_inputs.clear();
            delete setTunIn->m_children[0];
            setTunIn->m_children.erase(setTunIn->m_children.begin());
            tempVarPossTunIn->m_poss->DeleteNode(tempVarPossTunIn);
          }
          tun->m_pset->RemoveInTun(setTunIn);
          tun->RemoveAllInputs2Way();
          tun->m_poss->DeleteNode(tun);
          
          
          while (!setTunOut->m_inputs.empty()) {
            Node *tempVarPossTunOut = setTunOut->Input(0);
            tempVarPossTunOut->RemoveAllInputs2Way();
            delete tempVarPossTunOut->m_children[0];
            tempVarPossTunOut->m_children.clear();
            delete setTunOut->m_inputs[0];
            setTunOut->m_inputs.erase(setTunOut->m_inputs.begin());
            tempVarPossTunOut->m_poss->DeleteNode(tempVarPossTunOut);
          }
          tun->m_pset->RemoveOutTun(setTunOut);
          tun->m_poss->DeleteNode(setTunOut);
          return;
        }
        else
          throw;
      }
      else {
        if (!tun->m_pset)
          throw;
        if (tun->m_pset->IsShadow()) {
          throw;
        }
        if (tun->m_children.empty())
          throw;
        
        Tunnel *newSrcTunIn = NULL;
        // If the tempvar goes into this loop, find if the src for the tempvar
        // also goes into the loop in the same way (i.e., split the same way or
        // not) and use that. Otherwise, create a similar input to the loop.
        NodeConnVecIter iter = newSrc->m_children.begin();
        for(; !newSrcTunIn && iter != newSrc->m_children.end(); ++iter) {
          if ((*iter)->m_num == newSrcNum) {
            if ((*iter)->m_n->IsTunnel(SETTUNIN)) {
              Tunnel *inTun = (Tunnel*)((*iter)->m_n);
              if (inTun->m_pset == tun->m_pset) {
                if (inTun->IsLoopTunnel())
                  throw;
                newSrcTunIn = inTun;
              }
            }
          }
        }
        if (!newSrcTunIn) {
          //need to create in tuns
          throw;
        }
        //go through each poss tuns
        string name = tun->Child(0)->GetName(0).m_name;
        for (unsigned int i = 0; i < tun->m_children.size(); ++i) {
          Tunnel *tempVarPossTunIn = (Tunnel*)(tun->Child(i));
          Tunnel *newSrcPossTunIn = (Tunnel*)(newSrcTunIn->Child(i));
          if (!tempVarPossTunIn->m_children.empty()) {
            TempVarNode *newTemp = new TempVarNode(tempVarNode->m_info.GetDist(), name);
            newTemp->AddInput(newSrcPossTunIn, 0);
            if (newTemp->GetNameStr(0) != tempVarPossTunIn->GetNameStr(0)) {
              cout << newTemp->GetNameStr(0) << endl;
              cout << tempVarPossTunIn->GetNameStr(0) << endl;
              throw;
            }
            tempVarPossTunIn->RedirectChildren(0,newTemp,0);
	    if (newTemp->m_children.size() == 1) {
	      ((DLANode*)(newTemp->Child(0)))->AlignInfo(newTemp->m_align,
					   newTemp->m_alignModes,
					   newTemp->m_alignModesSrc);
	    }
            newSrcPossTunIn->m_poss->AddNode(newTemp);
	    newSrcPossTunIn->m_poss->m_fullyExpanded = false;
          }
        }
        
        while (!tun->m_children.empty()) {
          Node *tempVarPossTunIn = tun->Child(0);
          tempVarPossTunIn->RemoveAllChildren2Way();
          delete tempVarPossTunIn->m_inputs[0];
          tempVarPossTunIn->m_inputs.clear();
          delete tun->m_children[0];
          tun->m_children.erase(tun->m_children.begin());
          tempVarPossTunIn->m_poss->DeleteNode(tempVarPossTunIn);
        }
        tun->m_pset->RemoveInTun(tun);
        tun->RemoveAllInputs2Way();
        tun->m_poss->DeleteNode(tun);
        return;
      }
    }
    else
      throw;
  }
  else {
    throw;
  }
}

void MoveTempVarNodeIntoLoop::Apply(Node *node) const
{
  TempVarNode *tmp = (TempVarNode*)node;
  //MoveIn doesn't support this yet
  if (!tmp->m_sumDims.empty()) {
    throw;
  }
  if (!tmp->m_align.empty()) {
    throw;
  }
  MoveIn(tmp->Child(0), tmp->Input(0), tmp->InputConnNum(0), tmp);
  tmp->m_poss->DeleteChildAndCleanUp(tmp);
}

void TempVarNode::AddVariables(VarSet &set) const
{
  if (!m_align.empty()) {
    {
      Var var(ModeArrayVarType, m_alignModes);
      set.insert(var);
    }
    {
      Var var(ModeArrayVarType, m_alignModesSrc);
      set.insert(var);
    }
  }
}


void MoveTempVarNodeIntoSet::Apply(Node *node) const
{
  TempVarNode *tmp = (TempVarNode*)node;
  Poss *poss = tmp->m_poss;
  Tunnel *tun = (Tunnel*)(tmp->Child(0));
  RealPSet *set = (RealPSet*)(tun->m_pset);
  unsigned int oldTunNum = FindInTunVec(set->m_inTuns, tun);
  Node *input = tmp->Input(0);
  ConnNum num = tmp->InputConnNum(0);
  Tunnel *otherTun = NULL;
  for (auto child : input->m_children) {
    if (child->m_num == num && child->m_n->IsTunnel(SETTUNIN)) {
      if (((Tunnel*)(child->m_n))->m_pset == set) {
	otherTun = (Tunnel*)(child->m_n);
	break;
      }
    }
  }
  set->DisconnectFromSetsForMergingRecord();
  
  unsigned int tunNum;
  if (!otherTun) {
    otherTun = new Tunnel(SETTUNIN);
    otherTun->AddInput(input, num);
    otherTun->m_pset = set;
    poss->AddNode(otherTun);
    set->m_inTuns.push_back(otherTun);
    tunNum = set->m_inTuns.size() - 1;
    for(auto setPoss : set->m_posses) {
      Tunnel *possTun = new Tunnel(POSSTUNIN);
      possTun->AddInput(otherTun, 0);
      setPoss.second->m_inTuns.push_back(possTun);
      setPoss.second->AddNode(possTun);
      possTun->m_pset = set;
    }
  }
  else {
    tunNum = FindInTunVec(set->m_inTuns, otherTun);
  }

  for(auto setPoss : set->m_posses) {
    TempVarNode *newTemp = new TempVarNode;
    newTemp->Duplicate(tmp, true, false);
    setPoss.second->AddNode(newTemp);
    newTemp->AddInput(setPoss.second->m_inTuns[tunNum], 0);
    Node *oldPossTun = setPoss.second->m_inTuns[oldTunNum];
    oldPossTun->RedirectChildren(newTemp);
    setPoss.second->DeleteNode(oldPossTun);
  }
  
  set->m_inTuns.erase(set->m_inTuns.begin() + oldTunNum);
  poss->DeleteChildAndCleanUp(tun);
}

bool MoveTempVarNodeIntoSet::CanApply(const Node *node) const
{
  const TempVarNode *tmp = (TempVarNode*)node;
  if (tmp->m_children.size() == 1) {
    if (tmp->Child(0)->IsTunnel(SETTUNIN)) {
      const Tunnel *tun = (Tunnel*)(tmp->Child(0));
      BasePSet *set = tun->m_pset;
      if (set->IsLoop())
	return false;
      if (!set->IsReal())
	throw;
      if (!((RealPSet*)set)->m_shadows.empty())
	throw;
      return true;
    }
  }
  return false;
}


bool TempVarFromTempVar::CanApply(const Node *node) const
{
  const TempVarNode *tmp = (TempVarNode*)node;
  if (tmp->Input(0)->GetNodeClass() == TempVarNode::GetClass()) {
    const TempVarNode *inTmp = (TempVarNode*)(tmp->Input(0));
    if (!inTmp->m_sumDims.empty())
      throw;
    else
      return true;
  }
  else
    return false;
}

void TempVarFromTempVar::Apply(Node *node) const
{
  TempVarNode *tmp = (TempVarNode*)node;
  TempVarNode *inTmp = (TempVarNode*)(tmp->Input(0));
  
  tmp->ChangeInput2Way(inTmp, 0, 
		       inTmp->Input(0), inTmp->InputConnNum(0));
  
  if (tmp->m_align.empty() && !inTmp->m_align.empty()) {
    tmp->m_align = inTmp->m_align;
    tmp->m_alignModes = inTmp->m_alignModes;
    tmp->m_alignModesSrc = inTmp->m_alignModesSrc;
  }
  if (inTmp->m_children.empty()) {
    tmp->m_poss->DeleteChildAndCleanUp(inTmp);
  }
}



#endif //DOTENSORS
 
