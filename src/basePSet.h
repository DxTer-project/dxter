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



#pragma once

#include "base.h"
#include "poss.h"
#include "possTunnel.h"

class BasePSet
{
 public:
  NodeVec m_inTuns;
  NodeVec m_outTuns;
  Poss *m_ownerPoss;
  PossMMapIter m_currPoss;
  bool m_currHasPrinted;
  bool m_hasProped;
  bool m_isTopLevel;
  BasePSet();
  GraphNum NumPosses() {return m_posses.size();}
  virtual bool operator==(const PSet &rhs) const = 0;
  virtual void Prop() = 0;
  virtual void ClearBeforeProp();
  Node* InTun(unsigned int num) const;
  Node* OutTun(unsigned int num) const;
  bool GlobalSimplification(const TransMap &globalSimplifiers, const TransMap &simplifiers);
  virtual bool CanMerge(PSet *pset) const;
  virtual bool IsTransparent() const {return true;}
  bool MergePosses(const TransMap &simplifiers, CullFunction cullFunc);
  virtual GraphNum TotalCount() const = 0;
  virtual bool TakeIter(const TransMap &trans, const TransMap &simplifiers) = 0;
  virtual void InlinePoss(Poss *inliningPoss, PossMMap &newPosses) = 0;
  virtual void Duplicate(const BasePSet *orig, NodeMap &map, bool possMerging);
  virtual void ClearPrinted();
  virtual bool IsLoop() const {return false;}
  void RemoveInTun(Node *tun);
  void RemoveOutTun(Node *tun);
  virtual void FormSetAround() = 0;
  void ClearCurrPoss();
  bool IncrementCurrPoss();
  Cost EvalCurrPoss(TransConstVec &transList);
  Cost EvalAndSetBest();
  virtual void PrintCurrPoss(IndStream &out, GraphNum &graphNum);
  bool CanPrint() const;
  virtual Poss* GetCurrPoss() const = 0;
  void GetCurrTransVec(TransVec &transVec) const;
  void AddCurrPossVars(VarSet &set) const;

  void Flatten(ofstream &out) const;
  virtual void FlattenCore(ofstream &out) const {}
  void Unflatten(ifstream &in, SaveInfo &info);
  virtual void UnflattenCore(ifstream &in, SaveInfo &info) {}

  virtual void BuildDataTypeCache();
  virtual void ClearDataTypeCache();

#if DOBLIS
  virtual bool IsCritSect() const {return false;}
  Comm ParallelismWithinCurrentPosses() const;
  bool RemoveParallelization(Comm comm);
  void ReplaceAllComms(Comm comm1, Comm comm2);
#endif //DOBLIS

  const string& GetFunctionalityString() const;
};

