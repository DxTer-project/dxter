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



#pragma once

#define PRINTTRACKING 0
#define CHECKFORSETREUSE 0

#ifdef _OPENMP
#include "omp.h"
#endif
#include "base.h"
#include "poss.h"
#include <queue>
//#include "possTunnel.h"
#include "basePSet.h"

class Tunnel;
class ShadowPSet;


typedef map<int,int> IntMap;
typedef std::pair<int,int> IntMapPair;
typedef IntMap::iterator IntMapIter;
typedef IntMap::const_iterator IntMapConstIter;

/*
  this information tracks how two sets are fused w.r.t.
     connections between tunnels
  if there is an input of my set that is connected to the
     output of another set, the input tun number will map
     to the negative number of the tun of that set
  if there is an output for my set that is connected to 
     the input of another set, my negative output tun
     number is mapped to the positive in tun number of
     the other set
  if an input of mine is connected to the same node
     as the input of the other set, then the positive in tun
     number is mapped to the other set's positive in tun num
*/
class FusionInformation
{
 public:
  RealPSet *m_fused;
  IntMap m_map;
};


struct FusionInformationCompare {
  bool operator() (const FusionInformation& lhs, const FusionInformation& rhs) const{
    if (lhs.m_fused < rhs.m_fused)
      return true;
    else if (lhs.m_fused > rhs.m_fused)
      return false;
    else {
      if (lhs.m_map.size() < rhs.m_map.size())
	return true;
      else if (lhs.m_map.size() > rhs.m_map.size())
	return false;
      else {
	IntMapConstIter lhsIter = lhs.m_map.begin();
	IntMapConstIter rhsIter = rhs.m_map.begin();
	for(; lhsIter != lhs.m_map.end(); ++lhsIter, ++rhsIter) {
	  if (lhsIter->first < rhsIter->first) {
	    return true;
	  }
	  else if (lhsIter->first > rhsIter->first) {
	    return false;
	  }
	  else {
	    if (lhsIter->second < rhsIter->second)
	      return true;
	    else if (lhsIter->second > rhsIter->second)
	      return false;
	  }
	}
	return false;
      }
    }
  }
};

typedef std::pair<FusionInformation, RealPSet*> PSetMapPair;
typedef std::map<FusionInformation, RealPSet*, FusionInformationCompare> PSetMap;
typedef PSetMap::iterator PSetMapIter;
typedef PSetMap::const_iterator PSetMapConstIter;

class RealPSet : public BasePSet
{
 public:
#if DOTENSORS
  static RealPSetMMap m_setMap;
#ifdef _OPENMP
  static omp_lock_t m_lock;
#endif //_OPENMP
#endif
  PossMMap m_posses;
  string m_functionality;
  PSetVec m_shadows;
  RealPSet();
  RealPSet(Poss *poss);
  void Init(Poss *poss);
  PSetMap m_mergeMap;
  RealPSet *m_mergeLeft, *m_mergeRight;
  Cost m_cost;
  //after fusing loops, it's possible we have more input tuns connected to an input
  // than before since input of loop A could be output of loop B, but they split in different ways
  vector<vector<int>>  m_leftInMap, m_rightInMap;
  vector<int> m_leftOutMap, m_rightOutMap;
  virtual ~RealPSet();
  void UpdateRealPSetPointers(RealPSet *oldPtr, RealPSet *newPtr);
  void RemoveShadow(ShadowPSet *shadow);
  void AddPoss(Poss *poss);
  void AddPossesOrDispose(PossMMap &mmap, PossMMap *added = NULL);
  virtual GraphNum NumPosses() const {return m_posses.size();}
  bool operator==(const BasePSet &rhs) const;
  virtual Cost Prop();
  virtual bool TakeIter(const Universe *uni, int phase);
  virtual void ClearBeforeProp();
  virtual void Duplicate(const BasePSet *orig, NodeMap &map, bool possMerging, bool useShadows);
  void Migrate();
  void DisconnectFromSetsForMergingRecord();
  virtual BasePSet* GetNewInst() {return (BasePSet*)(new RealPSet);}
  virtual const PossMMap& GetPosses() const {return m_posses;}
  virtual PossMMap& GetPosses() {return m_posses;}
  virtual void PatchAfterDuplicate(NodeMap &map);
  void CombineAndRemoveTunnels();
  void RemoveAndDeletePoss(Poss *poss, bool removeFromMyList);
  void Simplify(const Universe *uni, int phase, bool recursive = false);
  void ClearFullyExpanded();
  virtual bool IsTransparent() const {return true;}
  void Cull(Phase phase);
  void Cull(CullFunction cullFunc);
  void CullWorstPerformers(double percentToCull, int ignoreThreshold);
  void CullAllBut(int num);
  bool MergePosses(const Universe *uni, int phase, CullFunction cullFunc);
  void FormSets(unsigned int phase);
  virtual GraphNum TotalCount() const;
  virtual void InlinePoss(Poss *inliningPoss, unsigned int num, PossMMap &newPosses);
  virtual ShadowPSet* GetNewShadow();
  virtual ShadowPSet* GetNewShadowDup(Poss *poss);

  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);

  virtual void BuildDataTypeCache();
  virtual void ClearDataTypeCache();

  void InlineAllSets();

#if DOTENSORS
  bool SamePSetWRTFunctionality(const RealPSet *other) const;
#endif

  bool RemoveLoops(bool *doneSomething);

#if DOBLIS
  Comm ParallelismWithinCurrentPosses() const;
  virtual bool RemoveParallelization(Comm comm);
#endif //DOBLIS

  virtual const string& GetFunctionalityString() const;
  virtual bool IsReal() const {return true;}
  virtual const RealPSet* GetReal() const {return this;}
  virtual RealPSet* GetReal() {return this;}
  static void GetFusionInformation(BasePSet *leftSet, BasePSet *rightSet,
				   RealPSet *realLeft, RealPSet *realRight,
				   FusionInformation &leftInfo, FusionInformation &rightInfo);
  static RealPSet* HasMergedWith(RealPSet *realLeft, RealPSet *realRight,
				 FusionInformation &leftInfo, FusionInformation &rightInfo);
  void SetDeletingRecursively();
  void ClearDeletingRecursively();

  bool EnforceMemConstraint(Cost costGoingIn, Cost maxMem, const StrSet &stillLive, Cost &highWater);
};


class PossCostComparison
{
 public:
  bool operator() (const Poss* lhs, const Poss* rhs) const
  {
    if (lhs->m_cost < 0)
      LOG_FAIL("replacement for throw call");
    if (rhs->m_cost < 0)
      LOG_FAIL("replacement for throw call");
    return lhs->m_cost > rhs->m_cost;
  }
};

typedef std::priority_queue<Poss*,std::vector<Poss*>,PossCostComparison> SortedPossQueue;
