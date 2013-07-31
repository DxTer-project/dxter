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



#include "loopSupport.h"
#include "distributions.h"
#include <cmath>
#include <climits>
#include "pack.h"
#include "critSect.h"

int Loop::M_currLabel = 0;

Size BSSizeToSize(BSSize size)
{
  switch(size)
  {
    case (USEELEMBS):
      return ELEM_BS;
    case (USEBLISMC):
      return BLIS_MC_BS;
    case (USEBLISKC):
      return BLIS_KC_BS;
    case (USEBLISNC):
      return BLIS_NC_BS;
    case (USEBLISOUTERBS):
      return BLIS_OUTER_BS;
    default:
      throw;
  }
}

string BSSizeToStr(BSSize size)
{
  switch(size)
  {
    case (USEELEMBS):
      throw;
    case (USEBLISMC):
      return "gemm_mc";
    case (USEBLISKC):
      return "gemm_kc";
    case (USEBLISNC):
      return "gemm_nc";
    case (USEBLISOUTERBS):
      return "bs_obj";
    default:
      throw;
  }
}

unsigned int GetNumElems(PartDir dir)
{
  switch(dir)
  {
    case(PARTRIGHT):
    case(PARTDOWN):
    case(PARTLEFT):
    case(PARTUPWARD):
      return 3;
    case (PARTDIAG):
    case (PARTDIAGBACK):
      return 9;
    default:
      cout << "Bad PartDir\n";
      return -1;
  }
}

string PartDirToStr(PartDir dir)
{
  switch(dir)
  {
    case(PARTRIGHT):
      return "right";
    case(PARTDOWN):
      return "down";
    case (PARTDIAG):
      return "diag";
    case(PARTLEFT):
      return "rightBack";
    case(PARTUPWARD):
      return "downBack";
    case (PARTDIAGBACK):
      return "diagBack";
    default:
      cout << "Bad PartDir\n";
      return "bad";
  }
}

Loop::Loop()
: PSet(), m_type(UNKNOWNLOOP), m_comm(CORECOMM), m_dim(BADDIM)
{
  AssignNewLabel();
  m_bsSize = BADBSSIZE;
}

Loop::Loop(LoopType type)
: m_type(type), m_comm(CORECOMM), m_dim(BADDIM)
{
  if (m_type == ELEMLOOP)
    m_bsSize = USEELEMBS;
  else
    m_bsSize = BADBSSIZE;
  AssignNewLabel();
}

Loop::Loop(LoopType type, Poss *poss, BSSize bsSize)
: PSet(poss), m_type(type), m_bsSize(bsSize), m_comm(CORECOMM), m_dim(BADDIM)
{
  unsigned int i;
  for(i = 0; i < poss->m_inTuns.size(); ++i) {
    if (poss->m_inTuns[i]->GetNodeClass() == Split::GetClass()) {
      if (((Split*)(poss->m_inTuns[i]))->m_isControlTun)
        ((Split*)(m_inTuns[i]))->m_isControlTun = true;
    }
  }
  AssignNewLabel();
  
  bool foundControl = false;
  NodeVecIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    Node *in = *iter;
    if (!in->IsLoopTunnel()) {
      cout << "non loop tunnel on loop!\n";
      throw;
    }
    if (in->GetNodeClass() == Split::GetClass())
    {
      Split *split = (Split*)in;
      
      if (split->m_isControlTun) {
        if (foundControl)
          throw;
        else
          foundControl = true;
      }
    }
  }
  if (!foundControl)
    throw;
}

void Loop::SanityCheck()
{
  PSet::SanityCheck();
  bool foundControl = false;
  NodeVecIter iter = m_inTuns.begin();
  //  cout << "****\n";
  for(; iter != m_inTuns.end(); ++iter) {
    Node *in = *iter;
    if (!in->IsLoopTunnel()) {
      cout << "non loop tunnel on loop!\n";
      throw;
    }
    if (in->GetNodeClass() == Split::GetClass()) {
      Split *split = (Split*)in;
      if (split->m_isControlTun) {
        if (foundControl)
          throw;
        else
          foundControl = true;
      }
    }
  }
  if (!foundControl)
    throw;
  iter = m_outTuns.begin();
  for(; iter != m_outTuns.end(); ++iter)
    if (!(*iter)->IsLoopTunnel()) {
      cout << "non loop tunnel on loop!\n";
      throw;
    }
}

bool Loop::CanMerge(PSet *pset) const
{
  //  return false;
  if (!pset->IsLoop())
    return false;
  if (m_bsSize != ((Loop*)pset)->m_bsSize)
    return false;
  Loop *loop = (Loop*)pset;
  if (loop->m_comm != CORECOMM || m_comm != CORECOMM)
    return false;
  if (m_type != loop->m_type)
    return false;
  const Split *split1 = GetControl();
  const Split *split2 = loop->GetControl();
  if (split1->NumberOfLoopExecs() != split2->NumberOfLoopExecs())
    return false;
  for(unsigned int i = 0; i < split1->NumberOfLoopExecs(); ++i) {
    if (split1->NumIters(i) != split2->NumIters(i))
      return false;
  }
  
#if DOSQOPHASE
  //If there's a PACKABLOCK in both sets, then we don't
  // want to fuse those loops.  We don't merge the packing
  // yet and they're likely using other packed B panels for the
  // computation.  That will make it impossible to separate the
  // different uses of the single packed B panel buffer
  bool found1 = false;
  bool found2 = false;
  for(unsigned int i = 0; !found1 && i < m_posses[0]->m_possNodes.size(); ++i) {
    Node *node = m_posses[0]->m_possNodes[i];
    if (node->GetNodeClass() == PackBuff::GetClass()) {
      PackBuff *buff = (PackBuff*)node;
      if (buff->m_packMat == PACKABLOCK)
        found1 = true;
    }
  }
  if (found1) {
    for(unsigned int i = 0; !found2 && i < loop->m_posses[0]->m_possNodes.size(); ++i) {
      Node *node = loop->m_posses[0]->m_possNodes[i];
      if (node->GetNodeClass() == PackBuff::GetClass()) {
        PackBuff *buff = (PackBuff*)node;
        if (buff->m_packMat == PACKABLOCK)
          found2 = true;
      }
    }
    if (found1 && found2)
      return false;
  }
#endif //DOSQOPHASE
  
  if (!PSet::CanMerge(pset))
    return false;
  
  const Loop *left=NULL, *right=NULL;
  //this is true if the left loop is actually on the left
  //otherwise, the order doesn't matter
  NodeVecConstIter iter = pset->m_inTuns.begin();
  for(; iter != pset->m_inTuns.end() && !left; ++iter) {
    const Node *inTun = *iter;
    for (unsigned int i = 0; i < inTun->m_inputs.size() && !left; ++i) {
      if (inTun->Input(i)->IsPossTunnel(SETTUNOUT) &&
          ((PossTunnel*)inTun->Input(i))->m_pset == this) {
        right = (Loop*)pset;
        left = this;
      }
    }
  }
  if (!left) {
    iter = m_inTuns.begin();
    for(; iter != m_inTuns.end() && !left; ++iter) {
      const Node *inTun = *iter;
      for (unsigned int i = 0; i < inTun->m_inputs.size() && !left; ++i) {
        if (inTun->Input(i)->IsPossTunnel(SETTUNOUT) &&
            ((PossTunnel*)inTun->Input(i))->m_pset == pset) {
          left = (Loop*)pset;
          right = this;
        }
      }
    }
  }
  
  if (left) {
    //make sure the output of the left loop is either only used as input to the
    // right loop or the right loop doesn't change it at all
    iter = left->m_outTuns.begin();
    for(; iter != left->m_outTuns.end(); ++iter) {
      LoopTunnel *leftOutTun = (LoopTunnel*)*iter;
      NodeConnVecConstIter connIter = leftOutTun->m_children.begin();
      bool foundConnToRightNotUpdated = false;
      
      bool foundConnNotToRight = false;
      for(; connIter != leftOutTun->m_children.end(); ++connIter) {
        const NodeConn *conn = *connIter;
        if (conn->m_n->IsLoopTunnel() && ((LoopTunnel*)conn->m_n)->m_pset == right) {
          LoopTunnel *rightInTun = (LoopTunnel*)conn->m_n;
          foundConnToRightNotUpdated |= !rightInTun->AllFullyUpdated();
          if (foundConnNotToRight) {
            if(foundConnToRightNotUpdated)
              return false;
          }
          //Check if the way the inputs/outputs are split are ok
          // for fusion
          if (leftOutTun->GetNodeClass() == Combine::GetClass()) {
            if (rightInTun->GetNodeClass() == Split::GetClass()) {
              if (((Combine*)leftOutTun)->m_dir != ((Split*)rightInTun)->m_dir) {
                if (!leftOutTun->IsConst() || !rightInTun->IsConst())
                  return false;
                return false;
                cout << "not yet supporting different directions";
                throw;
              }
            }
            else {
              if (!leftOutTun->IsConst() || !rightInTun->IsConst()) {
                return false;
              }
            }
          }
          else {
            if (rightInTun->GetNodeClass() == Split::GetClass()) {
              if (!leftOutTun->IsConst() || !rightInTun->IsConst()) {
                return false;
              }
            }
          }
        }
        else {
          foundConnNotToRight = true;
          if (foundConnToRightNotUpdated)
            return false;
        }
      }
    }
    
    //check we wouldn't have a cycle
    iter = right->m_outTuns.begin();
    for(; iter != right->m_outTuns.end(); ++iter) {
      const Node *outTun = *iter;
      for (unsigned int i = 0; i < outTun->m_children.size(); ++i) {
        if (outTun->Child(0)->IsPossTunnel(SETTUNIN) &&
            ((PossTunnel*)outTun->Child(0))->m_pset == left) {
          cout << "found cycle\n";
          throw;
        }
      }
    }
    
    
    //If the right loop uses a quadrant of a particular input
    // that is output by the left loop, make sure the left loop has
    // completely computed that quadrant
    iter = right->m_inTuns.begin();
    for(; iter != right->m_inTuns.end(); ++iter) {
      if (!(*iter)->IsLoopTunnel()) {
        throw;
      }
      const LoopTunnel *rightInTun = (LoopTunnel*)*iter;
      for (unsigned int i = 0; i < rightInTun->m_inputs.size() && left; ++i) {
        if (rightInTun->Input(0)->IsLoopTunnel() && rightInTun->Input(0)->IsPossTunnel(SETTUNOUT)) {
          const LoopTunnel *leftOutTun = (LoopTunnel*)(rightInTun->Input(0));
          if (((PossTunnel*)leftOutTun)->m_pset == left) {
            for(int quad = TL; quad < LASTQUAD; ++quad) {
              if (rightInTun->QuadInUse((Quad)quad, true)) {
                if (leftOutTun->GetUpStat((Quad)quad) != FULLUP) {
                  return false;
                }
              }
            }
          }
        }
      }
    }
    
    
    //If the left loop uses a quadrant of an output
    // used by the right loop, make sure the right loop
    // doesn't change it
    iter = left->m_inTuns.begin();
    for(; iter != left->m_inTuns.end(); ++iter) {
      if (!(*iter)->IsLoopTunnel()) {
        throw;
      }
      const LoopTunnel *leftInTun = (LoopTunnel*)*iter;
      const LoopTunnel *leftOutTun = leftInTun->GetMatchingOutTun();
      for(int quad = TL; quad < LASTQUAD; ++quad) {
        if (leftInTun->QuadInUse((Quad)quad,false)) {
          for (unsigned int i = 0; i < leftOutTun->m_children.size(); ++i) {
            if (leftOutTun->Child(i)->IsLoopTunnel()) {
              LoopTunnel *rightInTun = (LoopTunnel*)(leftOutTun->Child(i));
              if (rightInTun->m_pset == right) {
                if (!rightInTun->IsConst() && rightInTun->GetUpStat((Quad)quad) != NOTUP)
                  return false;
              }
            }
          }
        }
      }
    }
  }
  else {
    const PSet *leftSet = NULL;
    bool foundConnection = false;
    iter = pset->m_inTuns.begin();
    for (; iter != pset->m_inTuns.end(); ++iter) {
      LoopTunnel *tun1 = (LoopTunnel*)(*iter);
      Node *input = tun1->Input(0);
      unsigned int inNum = tun1->InputConnNum(0);
      NodeConnVecIter childIter = input->m_children.begin();
      for(; childIter != input->m_children.end(); ++childIter) {
        if ((*childIter)->m_num == inNum) {
          Node *child = (*childIter)->m_n;
          if (child->IsLoopTunnel() && ((LoopTunnel*)child)->m_pset == this) {
            LoopTunnel *tun2 = (LoopTunnel*)child;
            foundConnection = true;
            if (!tun1->IsConst() || !tun2->IsConst()) {
              
              
              if (tun1->IsConst()) {
                if (leftSet) {
                  if (leftSet != pset) {
                    cout << "!!!! Not supported\n";
                    return false;
                  }
                }
                else
                  leftSet = pset;
                LoopTunnel *outTun = tun1->GetMatchingOutTun();
                if (outTun) {
                  if (outTun->m_children.size()) {
                    cout << "out tun has children\n";
                    throw;
                  }
                }
              }
              else if (tun2->IsConst()){
                if (leftSet) {
                  if (leftSet != this) {
                    cout << "!!!! Not supported\n";
                    return false;
                  }
                }
                else
                  leftSet = this;
                LoopTunnel *outTun = tun2->GetMatchingOutTun();
                if (outTun) {
                  if (outTun->m_children.size()) {
                    cout << "out tun has children\n";
                    throw;
                  }
                }
              }
              else {
                cout << "not supported yet...\n";
                //Here, we might be able to fuse these two loops
                // that both use the same input...
                throw;
              }
              
              if (!leftSet)
                throw;
              
              LoopTunnel *leftTun;
              LoopTunnel *rightTun;
              
              if (leftSet == this) {
                leftTun = tun2;
                rightTun = tun1;
              }
              else {
                leftTun = tun1;
                rightTun = tun2;
              }
              
              for(int quad = TL; quad < LASTQUAD; ++quad) {
                if (leftTun->QuadInUse((Quad)quad,false)) {
                  if (!rightTun->IsConst() && rightTun->GetUpStat((Quad)quad) != NOTUP)
                    return false;
                }
              }
            }
          }
        }
      }
    }
    
    if (!foundConnection) {
      throw;
      return false;
    }
  }
  return true;
}

void Loop::PrintCurrPoss(IndStream &out, unsigned int &graphNum)
{
  string loopLevel = out.LoopLevel(1);
  
  NodeVecIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    if ((*iter)->GetNodeClass() == Split::GetClass()) {
      ((Split*)(*iter))->PrintVarDeclarations(out);
    }
  }
  if (m_type == BLISLOOP) {
    if (m_comm != CORECOMM)
      *out << "***Parallelized with communicator " << CommToStr(m_comm) << "; need correct output code\n";
    string idx = "idx" + loopLevel;
    string dimLen = "dimLen" + loopLevel;
    string bs = "bs" + loopLevel;
    
    Split *split = GetControl();
    
    string inputName = split->Input(0)->GetName(split->InputConnNum(0)).str();
    
    switch(split->m_dir) {
      case(PARTDIAGBACK):
      case(PARTDIAG):
        out.Indent();
        *out << dimLen << " = min( bli_obj_length_after_trans( " << inputName << " ), "
        << "bli_obj_width_after_trans( " << inputName << " ) );\n";
        break;
      case(PARTDOWN):
      case(PARTUPWARD):
        out.Indent();
        *out << dimLen << " = bli_obj_length_after_trans( " << inputName << " );\n";
        break;
      case (PARTRIGHT):
      case (PARTLEFT):
        out.Indent();
        *out << dimLen << " = bli_obj_width_after_trans( " << inputName << " );\n";
        break;
      default:
        throw;
    }
    out.Indent();
    *out << "for ( " << idx << " = 0; " << idx << " < " << dimLen << "; "
    << idx << " += " << bs <<" ) {\n";
    out.Indent(1);
    *out << bs;
    switch(split->m_dir) {
      case(PARTDOWN):
      case(PARTDIAG):
      case (PARTRIGHT):
        *out << " = bli_determine_blocksize_f( " ;
        break;
      case (PARTLEFT):
      case(PARTUPWARD):
      case(PARTDIAGBACK):
        *out << " = bli_determine_blocksize_b( " ;
        break;
      default:
        throw;
    }
    
    *out << idx << ", " << dimLen
    << ", &" << inputName << ", " << BSSizeToStr(m_bsSize) << " );\n";
    
    loopLevel = out.LoopLevel(2);
    idx = "idx" + loopLevel;
    dimLen = "dimLen" + loopLevel;
    bs = "bs" + loopLevel;
    out.Indent(1);
    *out << "dim_t " << idx << ", " << dimLen << ", " << bs << ";\n";
  }
  else {
    if (m_comm != CORECOMM)
      throw;
  }
  
  PSet::PrintCurrPoss(out, graphNum);
  if (m_type == BLISLOOP) {
    out.Indent();
    *out << "}\n";
  }
}

void Loop::Duplicate(const PSet *orig, NodeMap &map, bool possMerging)
{
  PSet::Duplicate(orig,map,possMerging);
  Loop *loop = (Loop*)orig;
  m_label = loop->m_label;
  m_bsSize = loop->m_bsSize;
  m_comm = loop->m_comm;
  m_dim = loop->m_dim;
  if (loop->m_bsSize >= BADBSSIZE) {
    cout << "duplicating a loop with zero blocksize\n";
    throw;
  }
}

void Loop::AssignNewLabel()
{
  m_label.insert(M_currLabel++);
}

bool Loop::WorthFusing(Loop *loop)
{
  NodeVecIter iter = m_outTuns.begin();
  for(; iter != m_outTuns.end(); ++iter) {
    Node *out = *iter;
    NodeConnVecIter iter2 = out->m_children.begin();
    for(; iter2 != out->m_children.end(); ++iter2) {
      Node *child = (*iter2)->m_n;
      if (child->IsPossTunnel(SETTUNIN)) {
        if (((PossTunnel*)child)->m_pset == loop)
          return true;
      }
    }
  }
  iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    Node *in = *iter;
    NodeConnVecIter iter2 = in->m_inputs.begin();
    for(; iter2 != in->m_inputs.end(); ++iter2) {
      Node *input = (*iter2)->m_n;
      if (input->IsPossTunnel(SETTUNOUT)) {
        if (((PossTunnel*)input)->m_pset == loop)
          return true;
      }
      unsigned int num = (*iter2)->m_num;
      NodeConnVecIter iter3 = input->m_children.begin();
      for(; iter3 != input->m_children.end(); ++iter3) {
        if ((*iter3)->m_num == num) {
          Node *otherChild = (*iter3)->m_n;
          if (otherChild->IsPossTunnel(SETTUNIN)) {
            if (((PossTunnel*)otherChild)->m_pset == loop)
              return true;
          }
        }
      }
    }
  }
  return false;
}

void Loop::SetBS(BSSize size)
{
  m_bsSize = size;
}

int Loop::GetBS() const
{
  return BSSizeToSize(m_bsSize);
}

Split* Loop::GetControl() const
{
  Split *control = NULL;
  NodeVecConstIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    Node *node = *iter;
    if (node->GetNodeClass() == Split::GetClass()) {
      Split *split = (Split*)node;
      if (split->m_isControlTun) {
        control = split;
      }
    }
  }
  if (!control)
    throw;
  return control;
}

void Loop::FlattenCore(ofstream &out) const
{
  WRITE(m_type);
  WRITE(m_bsSize);
  WRITE(m_comm);
  WRITE(m_dim);
  unsigned int size = m_label.size();
  WRITE(size);
  IntSetConstIter iter = m_label.begin();
  for(; iter != m_label.end(); ++iter)
    WRITE(*iter);
}


void Loop::UnflattenCore(ifstream &in, SaveInfo &info)
{
  READ(m_type);
  READ(m_bsSize);
  READ(m_comm);
  READ(m_dim);
  unsigned int size;
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
    int tmp;
    READ(tmp);
    m_label.insert(tmp);
  }
}

void Loop::FlattenStatic(ofstream &out)
{
  WRITE(M_currLabel);
}

void Loop::UnflattenStatic(ifstream &in)
{
  READ(M_currLabel);
}

void Loop::FillTunnelSizes()
{
  Split *control = GetControl();
  if (!control)
    throw;
  bool upToDate = true;
  NodeVecIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    LoopTunnel *tun = (LoopTunnel*)(*iter);
    if (!tun->m_msizes) {
      upToDate = false;
      break;
    }
  }
  if (upToDate)
    return;
  iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    LoopTunnel *tun = (LoopTunnel*)(*iter);
    tun->ClearSizeCache();
  }
  iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    Node *in = *iter;
    if (!in->IsLoopTunnel())
      throw;
    ((LoopTunnel*)in)->StartFillingSizes();
  }
  unsigned int numExecs = control->NumberOfLoopExecs();
  for(unsigned int i = 0; i < numExecs; ++i) {
    unsigned int numIters = control->NumIters(i);
    if (numIters) {
      iter = m_inTuns.begin();
      for (; iter != m_inTuns.end(); ++iter) {
        LoopTunnel *in = (LoopTunnel*)(*iter);
        in->AppendSizes(i, numIters, NumGroupsInComm(m_comm));
      }
    }
  }
}

void Loop::BuildSizeCache()
{
  FillTunnelSizes();
  PSet::BuildSizeCache();
}

LoopTunnel* Loop::CreateNewLoopTunnels(Node *input,
                                       unsigned int num, Poss *possToCareAbout, UpStat stat)
{
  LoopTunnel *newTunIn = new LoopTunnel(SETTUNIN);
  newTunIn->SetAllStats(stat);
  newTunIn->SetPSet(this);
  m_ownerPoss->AddNode(newTunIn);
  m_inTuns.push_back(newTunIn);
  
  LoopTunnel *newTunOut = new LoopTunnel(SETTUNOUT);
  newTunOut->SetAllStats(stat);
  newTunOut->SetPSet(this);
  m_ownerPoss->AddNode(newTunOut);
  m_outTuns.push_back(newTunOut);
  
  newTunIn->AddInput(input,num);
  
  Node *ret = NULL;
  PossVecIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter) {
    Poss *currPoss = *iter;
    LoopTunnel *newPossTunIn = new LoopTunnel(POSSTUNIN);
    newPossTunIn->SetAllStats(stat);
    currPoss->AddNode(newPossTunIn);
    currPoss->m_inTuns.push_back(newPossTunIn);
    newPossTunIn->AddInput(newTunIn);
    LoopTunnel *newPossTunOut = new LoopTunnel(POSSTUNOUT);
    newPossTunOut->SetAllStats(stat);
    currPoss->AddNode(newPossTunOut);
    currPoss->m_outTuns.push_back(newPossTunOut);
    newPossTunOut->AddInput(newPossTunIn,0);
    newPossTunOut->AddInput(newPossTunIn,1);
    newTunOut->AddInput(newPossTunOut,0);
    if(currPoss == possToCareAbout)
      ret = newPossTunIn;
  }
  if (!ret)
    throw;
  return (LoopTunnel*)ret;
}

void HardDeleteNode(Node *node)
{
  node->m_poss->RemoveFromGraphNodes(node);
  NodeConnVecIter iter = node->m_children.begin();
  for(; iter != node->m_children.end(); ++iter)
    delete *iter;
  node->m_children.clear();
  iter = node->m_inputs.begin();
  for(; iter != node->m_inputs.end(); ++iter)
    delete *iter;
  node->m_inputs.clear();
}

void Loop::TryToDeleteLoopTunnelSetAndCleanUp(LoopTunnel *tun)
{
  if (tun->GetNodeClass() != LoopTunnel::GetClass()) {
    cout << "only handling LoopTunnels specifically, not splits\n";
    throw;
  }
  if (tun->m_tunType != POSSTUNIN)
    throw;
  
  LoopTunnel *setTunIn = (LoopTunnel*)(tun->Input(0));
  
  //If there is an poss that uses the LoopTunnel's input,
  // then we shouldn't remove it from all
  //Check all children of the set input (poss inputs)
  // by checking all of their children - if all children's
  // children are PossTunOuts, then we're good.
  NodeConnVecIter childIter = setTunIn->m_children.begin();
  for(; childIter != setTunIn->m_children.end(); ++childIter) {
    Node *child = (*childIter)->m_n;
    NodeConnVecIter childIter2 = child->m_children.begin();
    for(; childIter2 != child->m_children.end(); ++childIter2) {
      Node *childChild = (*childIter2)->m_n;
      if (!childChild->IsPossTunnel(POSSTUNOUT))
        return;
    }
  }
  
  //At this point, we've checked, and there is no reason
  // to keep this set of loop tunnels around
  unsigned int tunNum = UINT_MAX;
  for(unsigned int i = 0; i < m_inTuns.size(); ++i) {
    if (m_inTuns[i] == setTunIn) {
      tunNum = i;
      break;
    }
  }
  if (tunNum == UINT_MAX)
    throw;
  
  Node *setTunOut = m_outTuns[tunNum];
  
  if (setTunOut->m_children.size() > 0)
    throw;
  
  NodeConnVecIter outIter = setTunOut->m_inputs.begin();
  for(; outIter != setTunOut->m_inputs.end(); ++outIter) {
    HardDeleteNode((*outIter)->m_n);
  }
  
  m_outTuns.erase(m_outTuns.begin()+tunNum);
  HardDeleteNode(setTunOut);
  
  childIter = setTunIn->m_children.begin();
  for(; childIter != setTunIn->m_children.end(); ++childIter) {
    Node *child = (*childIter)->m_n;
    HardDeleteNode(child);
    delete *childIter;
  }
  setTunIn->m_children.clear();
  
  m_inTuns.erase(m_inTuns.begin()+tunNum);
  
  setTunIn->m_poss->DeleteChildAndCleanUp(setTunIn, true);
}

void Loop::Parallelize(Comm comm)
{
  if (NumGroupsInComm(comm) <= 1)
    throw;
  if (RemoveParallelization(comm))
    throw;
  m_comm = comm;
  if (m_ownerPoss)
    m_ownerPoss->m_hasChecked = false;
  if (!HasIndepIters()) {
    bool found = false;
    //If this code changes, reflace in OnlyParallelizedOnNonIndependentData
    NodeVecIter iter = m_inTuns.begin();
    for(; iter != m_inTuns.end(); ++iter) {
      LoopTunnel *tun = (LoopTunnel*)(*iter);
      if (!tun->IndepIters() && !tun->InputIsTemp()) {
        if (found)
          throw;
        found = true;
        unsigned int numOut = tun->NumOutputs();
        if (tun->GetNodeClass() == Split::GetClass())
          --numOut;
        int i;
        for(i = 0; i < (int)(tun->m_children.size()); ++i) {
          Node *possTun = (tun->m_children[i])->m_n;
          Poss *poss = possTun->m_poss;
          NodeSet nodeSet;
          for (unsigned int j = 0; j < numOut; ++j) {
            AddUsersOfLiveOutput(possTun, j, nodeSet);
          }
          if (!nodeSet.size())
            throw;
          poss->FillClique(nodeSet);
	  poss->m_hasChecked = false;
          CritSect *crit = (CritSect*)(poss->FormSetForClique(nodeSet, true));
          if (crit->RemoveParallelization(CORECOMM)) {
            //This critical section is around some hierarchy of PSets
            // from which parallel code cannot be removed without getting
            // rid of all code
            RemoveAndDeletePoss(poss, true);
            --i;
          }
        }
      }
    }
  }

  
  
  ClearSizeCache();
  
  //If we're parallelizing a loop that is on a poss
  // that just got duplicated as part of a transformation,
  // then that duplicated poss doesn't have its size cache.
  //We need to form it and this loop's size cache (which will be
  // different thanks to paralellization.
  if (m_ownerPoss)
    m_ownerPoss->BuildSizeCache();
}

bool ContainsOnlyParallelization(PSet *set)
{
  if (set->IsLoop()) {
    Loop *loop = (Loop*)set;
    if (loop->IsParallel()) {
      return true;
    }
  }
  
  PossVecIter iter = set->m_posses.begin();
  for(; iter != set->m_posses.end(); ++iter) {
    Poss *poss = *iter;
    bool foundParSet = false;
    PSetVecIter setIter = poss->m_sets.begin();
    for(; !foundParSet && setIter != poss->m_sets.end(); ++setIter) {
      PSet *set = *setIter;
      if (ContainsOnlyParallelization(set))
        foundParSet = true;
    }
    if (!foundParSet) {
      bool foundParNode = false;
      NodeVecIter nodeIter = poss->m_possNodes.begin();
      for(; !foundParNode && nodeIter != poss->m_possNodes.end(); ++nodeIter) {
        if ((*nodeIter)->IsParallel()) {
          foundParNode = true;
        }
      }
      if (!foundParNode)
        return false;
    }
  }
  return true;
}

bool ContainsOnlyParallelization(const NodeSet &set)
{
  PSetSet setSet;
  NodeSetIter iter = set.begin();
  for(; iter != set.end(); ++iter) {
    Node *node = *iter;
    if (node->IsParallel())
      return true;
    if (node->IsPossTunnel(SETTUNIN)) {
      PossTunnel *tun = (PossTunnel*)node;
      if (setSet.find(tun->m_pset) == setSet.end()) {
        setSet.insert(tun->m_pset);
        if (ContainsOnlyParallelization(tun->m_pset))
          return true;
      }
    }
  }
  return false;
}

bool Loop::OnlyParallelizedOnNonIndependentData() const
{
  NodeVecConstIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    //If this code changes, reflect in Parallelize
    const LoopTunnel *tun = (LoopTunnel*)(*iter);
    if (!tun->IndepIters() && !tun->InputIsTemp()) {
      unsigned int numOut = tun->NumOutputs();
      if (tun->GetNodeClass() == Split::GetClass())
        --numOut;
      bool foundNonPar = false;
      NodeConnVecConstIter connIter = tun->m_children.begin();
      for( ; !foundNonPar && connIter != tun->m_children.end(); ++ connIter) {
        Node *possTun = (*connIter)->m_n;
        Poss *poss = possTun->m_poss;
        NodeSet nodeSet;
        for (unsigned int i = 0; i < numOut; ++i) {
          AddUsersOfLiveOutput(possTun, i, nodeSet);
        }
        if (!nodeSet.size())
          throw;
        poss->FillClique(nodeSet);
        if (!ContainsOnlyParallelization(nodeSet)) {
          foundNonPar = true;
        }
      }
      if (!foundNonPar)
        return true;
    }
  }
  return false;
}

bool Loop::HasIndepIters() const
{
  NodeVecConstIter iter = m_inTuns.begin();
  for (; iter != m_inTuns.end(); ++iter) {
    const LoopTunnel *in = (LoopTunnel*)(*iter);
    if (!in->IndepIters()) {
      return false;
    }
  }
  return true;
}
