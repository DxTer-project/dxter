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



#pragma once

#include "base.h"
#include "poss.h"
#include "possTunnel.h"

class PSet
{
 public:
  PossMMap m_posses;
  NodeVec m_inTuns;
  NodeVec m_outTuns;
  bool m_hasProped;
  bool m_isTopLevel;
  Poss *m_ownerPoss;
  PossMMapIter m_currPoss;
  PSet();
  PSet(Poss *poss);
  virtual ~PSet();
  void AddPoss(Poss *poss);
  void AddPossesOrDispose(PossMMap &mmap, PossMMap *added = NULL);
  unsigned int NumPosses() {return m_posses.size();}
  virtual void SanityCheck();
  bool operator==(const Poss &rhs) const;
  bool operator==(const PSet &rhs) const;
  virtual void Prop();
  void Cull(Phase phase);
  Node* InTun(unsigned int num) const;
  Node* OutTun(unsigned int num) const;
  void ClearBeforeProp();
  bool TakeIter(const TransMap &trans, const TransMap &simplifiers);
  bool GlobalSimplification(const TransMap &globalSimplifiers, const TransMap &simplifiers);
  unsigned int Size() { return m_posses.size();} 
  virtual void Duplicate(const PSet *orig, NodeMap &map, bool possMerging);
  virtual PSet* GetNewInst() {return new PSet;}
  void PatchAfterDuplicate(NodeMap &map);
  void CombineAndRemoveTunnels();
  void RemoveAndDeletePoss(Poss *poss, bool removeFromMyList);
  void Simplify(const TransMap &simplifiers);
  //  void RemoveDups();
  void ClearFullyExpanded();
  virtual bool CanMerge(PSet *pset) const;
  virtual bool IsTransparent() const {return true;}
  bool MergePosses(const TransMap &simplifiers, CullFunction cullFunc);
  void FormSets(unsigned int phase);
  unsigned int TotalCount() const;
  void InlinePoss(Poss *inliningPoss, PossMMap &newPosses);
  virtual void ClearPrinted();
  virtual bool IsLoop() const {return false;}
  virtual bool IsCritSect() const {return false;}
  void RemoveInTun(Node *tun);
  void RemoveOutTun(Node *tun);
  void Cull(CullFunction cullFunc);
  void FormSetAround();
  void ClearCurrPoss();
  bool IncrementCurrPoss();
  Cost EvalCurrPoss(TransConstVec &transList);
  virtual void PrintCurrPoss(IndStream &out, unsigned int &graphNum);
  bool CanPrint() const;
  Poss* GetCurrPoss() const;
  void GetCurrTransVec(TransVec &transVec) const;
  Comm ParallelismWithinCurrentPosses() const;

  void Flatten(ofstream &out) const;
  virtual void FlattenCore(ofstream &out) const {}
  void Unflatten(ifstream &in, SaveInfo &info);
  virtual void UnflattenCore(ifstream &in, SaveInfo &info) {}

  virtual void BuildSizeCache();
  virtual void ClearSizeCache();
  bool RemoveParallelization(Comm comm);
  void ReplaceAllComms(Comm comm1, Comm comm2);
};

