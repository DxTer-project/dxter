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

#define PRINTTRACKING 0

#include "base.h"
#include "poss.h"
#include <queue>
//#include "possTunnel.h"
#include "basePSet.h"

class Tunnel;
class ShadowPSet;

class RealPSet : public BasePSet
{
 public:
  PossMMap m_posses;
  string m_functionality;
  PSetVec m_shadows;
  RealPSet();
  RealPSet(Poss *poss);
  void Init(Poss *poss);
  PSetMap m_mergeMap;
  RealPSet *m_mergeLeft, *m_mergeRight;
  Cost m_cost;
  vector<int>  m_leftInMap, m_rightInMap, m_leftOutMap, m_rightOutMap;
  virtual ~RealPSet();
  void UpdateRealPSetPointers(RealPSet *oldPtr, RealPSet *newPtr);
  void RemoveShadow(ShadowPSet *shadow);
  void AddPoss(Poss *poss);
  void AddPossesOrDispose(PossMMap &mmap, PossMMap *added = NULL);
  virtual GraphNum NumPosses() const {return m_posses.size();}
  bool operator==(const BasePSet &rhs) const;
  virtual Cost Prop();
  virtual bool TakeIter(const TransMap &trans, const TransMap &simplifiers);
  virtual void ClearBeforeProp();
  virtual void Duplicate(const BasePSet *orig, NodeMap &map, bool possMerging, bool useShadows);
  void Migrate();
  virtual BasePSet* GetNewInst() {return (BasePSet*)(new RealPSet);}
  virtual const PossMMap& GetPosses() const {return m_posses;}
  virtual PossMMap& GetPosses() {return m_posses;}
  virtual void PatchAfterDuplicate(NodeMap &map);
  void CombineAndRemoveTunnels();
  void RemoveAndDeletePoss(Poss *poss, bool removeFromMyList);
  void Simplify(const TransMap &simplifiers, bool recursive = false);
  void ClearFullyExpanded();
  virtual bool IsTransparent() const {return true;}
  void Cull(Phase phase);
  void Cull(CullFunction cullFunc);
  void CullWorstPerformers(double percentToCull, int ignoreThreshold);
  bool MergePosses(const TransMap &simplifiers, CullFunction cullFunc);
  void FormSets(unsigned int phase);
  virtual GraphNum TotalCount() const;
  virtual void InlinePoss(Poss *inliningPoss, PossMMap &newPosses);
  virtual ShadowPSet* GetNewShadow();
  virtual ShadowPSet* GetNewShadowDup(Poss *poss);

  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);

  virtual void BuildDataTypeCache();
  virtual void ClearDataTypeCache();

  void InlineAllSets();

  void SetInTunsAsPrinted();
  bool RemoveLoops(bool *doneSomething);

#if DOBLIS
  Comm ParallelismWithinCurrentPosses() const;
  virtual bool RemoveParallelization(Comm comm);
#endif //DOBLIS

  virtual const string& GetFunctionalityString() const;
  virtual bool IsReal() const {return true;}
  virtual const RealPSet* GetReal() const {return this;}
  virtual RealPSet* GetReal() {return this;}

  RealPSet* HasMergedWith(RealPSet *set, bool checkOtherOrder=true);
  void SetDeletingRecursively();
  void ClearDeletingRecursively();
};


class PossCostComparison
{
 public:
  bool operator() (const Poss* lhs, const Poss* rhs) const
  {
    if (lhs->m_cost < 0)
      throw;
    if (rhs->m_cost < 0)
      throw;
    return lhs->m_cost > rhs->m_cost;
  }
};

typedef std::priority_queue<Poss*,std::vector<Poss*>,PossCostComparison> SortedPossQueue;
