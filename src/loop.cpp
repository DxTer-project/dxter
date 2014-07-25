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



#include "loopSupport.h"
#include "elemRedist.h"
#include <cmath>
#include <climits>
#include "pack.h"
#include "critSect.h"
#include "blis.h"

int Loop::M_currLabel = 0;

BSSize LLDLAMu(USELLDLAMU);
BSSize LLDLA2Mu(USELLDLA2MU);
BSSize LLDLA3Mu(USELLDLA3MU);
BSSize BadBS(BADBSSIZE);
BSSize UnitBS(USEUNITBS);

string BSSizeToVarName(BSSize size)
{
  switch(size.m_val)
    {
#if DOLLDLA
    case (USELLDLAMU):
      return MU_VAR_NAME;
    case (USELLDLA2MU):
      return "(2*" + (string)(MU_VAR_NAME) + ")";
    case (USELLDLA3MU):
      return "(3*" + (string)(MU_VAR_NAME) + ")";
#endif
    default:
      throw;
    }
}



string BSSizeToStr(BSSize size)
{
  switch(size.m_val)
  {
#if DOELEM
    case (USEELEMBS):
      throw;
#elif DOBLIS
    case (USEBLISMC):
      return "gemm_mc";
    case (USEBLISKC):
      return "gemm_kc";
    case (USEBLISNC):
      return "gemm_nc";
    case (USEBLISOUTERBS):
      return "bs_obj";
#endif
    default:
      throw;
  }
}

string BSSizeToSubSizeStr(BSSize size)
{
  switch(size.m_val)
  {
#if DOELEM
    case (USEELEMBS):
      throw;
#elif DOBLIS
    case (USEBLISMC):
      return "gemm_mr";
    case (USEBLISKC):
      return "gemm_kr";
    case (USEBLISNC):
      return "gemm_nr";
    case (USEBLISOUTERBS):
      throw;
#endif
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
: PSet(), 
#if TWOD
 m_dim(BADDIM),
#endif
  m_type(UNKNOWNLOOP)
#if DOBLIS
, m_comm(CORECOMM)
#endif
{
  AssignNewLabel();
  m_bsSize = BadBS;
}

Loop::Loop(LoopType type)
: 
#if TWOD
  m_dim(BADDIM),
#endif
  m_type(type)
#if DOBLIS
, m_comm(CORECOMM)
#endif

{
#if DOELEM
  if (m_type == ELEMLOOP)
    m_bsSize = USEELEMBS;
  else
#endif
    m_bsSize = BadBS;
  AssignNewLabel();
}

Loop::Loop(LoopType type, Poss *poss, BSSize bsSize)
: PSet(poss)
#if TWOD
, m_dim(BADDIM)
#endif
, m_type(type), m_bsSize(bsSize)
#if DOBLIS
, m_comm(CORECOMM)
#endif

{
  unsigned int i;
  for(i = 0; i < poss->m_inTuns.size(); ++i) {
    LoopTunnel *tun = (LoopTunnel*)(poss->m_inTuns[i]);
    if (tun->IsSplit()) {
      SplitBase *split = (SplitBase*)tun;
      if (split->m_isControlTun)
        split->m_isControlTun = true;
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
    if (((LoopTunnel*)in)->IsSplit())
    {
      SplitBase *split = (SplitBase*)in;
      
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

void Loop::Prop()
{
  PSet::Prop();

  
  bool foundControl = false;
  NodeVecIter iter = m_inTuns.begin();
  //  cout << "****\n";
  for(; iter != m_inTuns.end(); ++iter) {
    Node *in = *iter;
    if (!in->IsLoopTunnel()) {
      cout << "non loop tunnel on loop!\n";
      throw;
    }
    if (((LoopTunnel*)in)->IsSplit()) {
      SplitBase *split = (SplitBase*)in;
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
  if (!pset->IsLoop())
    return false;
  if (m_bsSize != ((Loop*)pset)->m_bsSize)
    return false;
  Loop *loop = (Loop*)pset;
#if DOBLIS
  if (loop->m_comm != CORECOMM || m_comm != CORECOMM)
    return false;
#endif
  if (m_type != loop->m_type)
    return false;
  const SplitBase *splitBase1 = GetControl();
  const SplitBase *splitBase2 = loop->GetControl();
  if (splitBase1->GetNodeClass() != SplitSingleIter::GetClass() ||
      splitBase2->GetNodeClass() != SplitSingleIter::GetClass())
    return false;
  const SplitSingleIter *split1 = (SplitSingleIter*)splitBase1;
  const SplitSingleIter *split2 = (SplitSingleIter*)splitBase2;
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
    for (ConnNum i = 0; i < inTun->m_inputs.size() && !left; ++i) {
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
      for (ConnNum i = 0; i < inTun->m_inputs.size() && !left; ++i) {
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
          if (leftOutTun->GetNodeClass() == CombineSingleIter::GetClass()) {
            if (rightInTun->GetNodeClass() == SplitSingleIter::GetClass()) {
#if TWOD
              if (((CombineSingleIter*)leftOutTun)->m_dir != ((SplitSingleIter*)rightInTun)->m_dir) {
#else
              if (((CombineSingleIter*)leftOutTun)->m_partDim != ((SplitSingleIter*)rightInTun)->m_partDim) {
#endif
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
            if (rightInTun->GetNodeClass() == SplitSingleIter::GetClass()) {
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
      for (ConnNum i = 0; i < rightInTun->m_inputs.size() && left; ++i) {
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
      ConnNum inNum = tun1->InputConnNum(0);
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

void Loop::PrintCurrPoss(IndStream &out, GraphNum &graphNum)
{
  Poss *poss = GetCurrPoss();
  NodeVecIter iter = poss->m_inTuns.begin();
  for(; iter != poss->m_inTuns.end(); ++iter) {
    if (((LoopTunnel*)(*iter))->IsSplit()) {
      ((SplitBase*)(*iter))->PrintVarDeclarations(out);
    }
  }

#if DOBLIS
  string loopLevel = out.LoopLevel(1);
  
  if (m_type == BLISLOOP) {
    if (m_comm != CORECOMM) {
      bool barrier = false;
      iter = m_inTuns.begin();
      for(; iter != m_inTuns.end() && !barrier; ++iter) {
	if (!FoundBarrier(*iter, 0, m_comm))
	  barrier = true;
      }
      if (barrier) {
	out.Indent();
	*out << "th_barrier( " << CommToStr(m_comm) << " );"
	     << "\t //barrier for dependency\n";

      }
      out.Indent();
      *out << "//// ***Parallelized with communicator "
	   << CommToStr(m_comm) << endl;
    }
    string idx = "idx" + loopLevel;
    string dimLen = "dimLen" + loopLevel;
    string bs = "bs" + loopLevel;
    
    SplitBase *splitBase = GetControl();
    if (splitBase->GetNodeClass() != SplitSingleIter::GetClass())
      throw;
    SplitSingleIter *split = (SplitSingleIter*)splitBase;
    
    string inputName = split->Input(0)->GetName(split->InputConnNum(0)).str();
    
    if (loopLevel == "1") {
      out.Indent();
      *out << "dim_t " << idx << ", " << dimLen << ", " << bs << ";\n";
    }
    
    
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
    if (m_comm != CORECOMM) {
      *out << idx << " = 0;\n";
      out.Indent();
      *out << "th_shift_start_end(&" << idx << ", &" << dimLen << ", "
      << CommToStr(GetSubComm(m_comm)) << ", "
      << "bli_blksz_for_obj( &" << split->GetInputNameStr(0)
      << ", " << BSSizeToSubSizeStr(m_bsSize) << "));\n";
      out.Indent();
      *out << "for ( ; " << idx << " < " << dimLen << "; "
      << idx << " += " << bs <<" ) {\n";
    }
    else {
#if DOSM
      Comm outerComm = split->WithinParallelism();
      Comm innerComm = ParallelismWithinCurrentPosses();
      if ((innerComm != ALLPROCCOMM && innerComm != ALLL2COMM) &&
	  (outerComm == CORECOMM || innerComm != GetSubComm(outerComm))) {
        if (innerComm == CORECOMM) {
          if (outerComm == CORECOMM)
            *out << "if (th_global_thread_id() != 0)\n";
          else
            *out << "if (th_thread_id( " << CommToStr(GetSubComm(outerComm)) << " ) != 0)\n";
        }
        else if (outerComm != CORECOMM && innerComm != GetSubComm(GetSubComm(outerComm))) {
          outerComm = split->WithinParallelism();
          innerComm = ParallelismWithinCurrentPosses();
          throw;
        }
        else {  
          *out << "if (th_group_id( " << CommToStr(innerComm) << " ) != 0)\n";
        }
        out.Indent(1);
        *out << dimLen << " = 0;\n";
        out.Indent();
      }
#endif // DOSM
      *out << "for ( " << idx << " = 0; " << idx << " < " << dimLen << "; "
      << idx << " += " << bs <<" ) {\n";
    }
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
#elif DOLLDLA
  if (m_type != LLDLALOOP)
    throw;
  if (m_bsSize != LLDLAMu &&
      m_bsSize != LLDLA2Mu &&
      m_bsSize != LLDLA3Mu)
    throw;
  SplitBase *split = GetControl();
  switch(m_dim) 
    {
    case (DIMM):
      {
	out.Indent();
	*out << "//Dim-m loop\n";;
	break;
      }
    case (DIMN):
      {
	out.Indent();
	*out << "//Dim-n loop\n";
	break;
      }
    case (DIMK):
      {
	out.Indent();
	*out << "//Dim-k loop\n";
	break;
    }
    default:
      break;
    }

  string loopLevel = split->GetLoopLevel();
  string lcv = "lcv" + loopLevel;
  out.Indent();
  *out << "for( " << lcv << " = ";
  
  if (split->GetNodeClass() != SplitSingleIter::GetClass())
    throw;

  bool needMin = false;
  if (split->m_dir == PARTDOWN) {
    *out << split->InputDataType(0).m_numRowsVar;
    if (!split->GetInputM(0)->EvenlyDivisibleBy(m_bsSize.Size()))
      needMin = true;
  }
  else if (split->m_dir == PARTRIGHT) {
    *out << split->InputDataType(0).m_numColsVar;
    if (!split->GetInputN(0)->EvenlyDivisibleBy(m_bsSize.Size()))
      needMin = true;
  }
  else
    throw;
  
  *out << "; " << lcv << " > 0; " << lcv << " -= ";

  
  *out << BSSizeToVarName(m_bsSize) << " ) {\n";

  out.Indent(1);
  *out << "const unsigned int num";
  if (split->m_dir == PARTDOWN) {
    *out << "Rows";
  }
  else if (split->m_dir == PARTRIGHT) {
    *out << "Cols";
  }
  else
    throw;

  if (needMin)
    *out << loopLevel << " = min( " << lcv << ", ";
  else
    *out << loopLevel << " = ( ";
  *out << BSSizeToVarName(m_bsSize) << " );\n";
#endif
  
  PSet::PrintCurrPoss(out, graphNum);

  iter = poss->m_inTuns.begin();
  for(; iter != poss->m_inTuns.end(); ++iter) {
    LoopTunnel *tun = (LoopTunnel*)(*iter);
    if (tun->IsSplit()) {
      SplitBase *split = (SplitBase*)tun;
      split->PrintIncrementAtEndOfLoop(m_bsSize, out);
    }
  }

#if DOBLIS||DOLLDLA
  if (m_type == BLISLOOP || m_type == LLDLALOOP) {
    out.Indent();
    *out << "}\n";
  }
#endif //DOBLIS||DOLLDLA
}
 
unsigned int Loop::LoopLevel() const
{
  unsigned int level = 0;
  Poss *poss = m_ownerPoss;
  while (poss) {
    if (!poss->m_pset)
      throw;
    if (poss->m_pset->IsLoop())
      ++level;
    poss = poss->m_pset->m_ownerPoss;
  }
  return level;
}
 
void Loop::Duplicate(const PSet *orig, NodeMap &map, bool possMerging)
{
  PSet::Duplicate(orig,map,possMerging);
  Loop *loop = (Loop*)orig;
  m_label = loop->m_label;
  m_bsSize = loop->m_bsSize;
#if DOBLIS
  m_comm = loop->m_comm;
#endif
#if TWOD
  m_dim = loop->m_dim;
#endif
}

void Loop::AssignNewLabel()
{
  m_label.insert(M_currLabel++);
}

bool Loop::WorthFusing(Loop *loop)
{
  //If we fuse two inner-most BLIS loops, there
  // could be two packed BPANEL buffers input into the same
  // loop, which means the buffers cannot be named the same.
  //Therefore, prohibit such fusion to allow for the same
  // buffer to be reused.
#if DOBLIS
  if (FindOtherPackBuffs((*(m_posses.begin())).second, PACKABLOCK, NULL)) {
    return false;
  }
#endif
  
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
      ConnNum num = (*iter2)->m_num;
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
  return (int)(m_bsSize.Size());
}

SplitBase* Loop::GetControl() const
{
  SplitBase *control = NULL;
  NodeVecConstIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    LoopTunnel *node = (LoopTunnel*)(*iter);
    if (node->IsSplit()) {
      SplitBase *split = (SplitBase*)node;
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
#if DOBLIS
  WRITE(m_comm);
#endif
#if TWOD
  WRITE(m_dim);
#endif
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
#if DOBLIS
  READ(m_comm);
#endif
#if TWOD
  READ(m_dim);
#endif
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
  SplitBase *control = GetControl();
  if (!control)
    throw;
  bool upToDate = true;
  NodeVecIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    LoopTunnel *tun = (LoopTunnel*)(*iter);
#if TWOD
    if (!tun->m_msizes) {
#else
    if (!tun->m_sizes) {
#endif
      upToDate = false;
      break;
    }
  }
  if (upToDate)
    return;
  iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    LoopTunnel *tun = (LoopTunnel*)(*iter);
    tun->ClearDataTypeCache();
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
#if DOBLIS
        in->AppendSizes(i, numIters, NumGroupsInComm(m_comm));
#else
        in->AppendSizes(i, numIters, 1);
#endif
      }
    }
  }
#if DODM
  iter = m_inTuns.begin();
  for (; iter != m_inTuns.end(); ++iter) {
    LoopTunnel *in = (LoopTunnel*)(*iter);
    in->UpdateLocalSizes();
  }
#endif
}

void Loop::BuildDataTypeCache()
{
  FillTunnelSizes();
  PSet::BuildDataTypeCache();
}

LoopTunnel* Loop::CreateNewLoopTunnels(Node *input,
                                       ConnNum num, Poss *possToCareAbout, UpStat stat)
{
  /*
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
    Poss *currPoss = (*iter).second;
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
  */
  throw;
  return NULL;
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

#if DOBLIS
void Loop::Parallelize(Comm comm)
{
  if (NumGroupsInComm(comm) <= 1)
    throw;
  if (RemoveParallelization(comm))
    throw;
  m_comm = comm;
  m_hasProped = false;
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
        if (tun->IsSplit()) {
	  if (tun->GetNodeClass() != SplitSingleIter::GetClass())
	    throw;
	  --numOut;
	}
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
	  throw;
#if 0
          CritSect *crit = (CritSect*)(poss->FormSetForClique(nodeSet, true));
          if (crit->RemoveParallelization(CORECOMM)) {
            //This critical section is around some hierarchy of PSets
            // from which parallel code cannot be removed without getting
            // rid of all code
            RemoveAndDeletePoss(poss, true);
            --i;
          }
          if (HasParallelCode(crit->m_posses[0])) {
            throw;
          }
#endif
        }
      }
    }
  }
  
  
  
  ClearDataTypeCache();
  
  //If we're parallelizing a loop that is on a poss
  // that just got duplicated as part of a transformation,
  // then that duplicated poss doesn't have its size cache.
  //We need to form it and this loop's size cache (which will be
  // different thanks to paralellization.
  if (m_ownerPoss)
    m_ownerPoss->BuildDataTypeCache();
}

bool ContainsOnlyParallelization(PSet *set)
{
  if (set->IsLoop()) {
    Loop *loop = (Loop*)set;
    if (loop->IsParallel()) {
      return true;
    }
  }
  
  PossMMapIter iter = set->m_posses.begin();
  for(; iter != set->m_posses.end(); ++iter) {
    Poss *poss = (*iter).second;
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
      if (tun->IsSplit()) {
	if (tun->GetNodeClass() != SplitSingleIter::GetClass())
	  throw;
	--numOut;
      }
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
#endif

#if TWOD
 void Loop::SetDimName(DimName dim)
 {
   m_functionality += (char)(48+dim);
   m_dim = dim;
 }
#endif



