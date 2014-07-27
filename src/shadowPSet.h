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
#include "tunnel.h"
#include "basePSet.h"

class ShadowPSet : public BasePSet
{
 public:
  RealPSet *m_realPSet;
  ShadowPSet();
  ShadowPSet(Poss *poss) {throw;}
  virtual ~ShadowPSet();
  virtual GraphNum NumPosses() const {return m_realPSet->NumPosses();}
  bool operator==(const BasePSet &rhs) const;
  virtual void Prop();
  virtual BasePSet* GetNewInst() {return new ShadowPSet;}
  void RemoveAndDeletePoss(Poss *poss, bool removeFromMyList);
  virtual bool IsTransparent() const {return true;}
  virtual GraphNum TotalCount() const;
  virtual const PossMMap& GetPosses() const {return m_realPSet->m_posses;}
  virtual PossMMap& GetPosses() {return m_realPSet->m_posses;}
  virtual void PatchAfterDuplicate(NodeMap &map) {}
  virtual void Duplicate(const BasePSet *orig, NodeMap &map, bool possMerging, bool useShadows);

  //  virtual void FlattenCore(ofstream &out) const {}
  //  virtual void UnflattenCore(ifstream &in, SaveInfo &info) {}

  virtual void BuildDataTypeCache();
  virtual void ClearDataTypeCache();

#if DOBLIS
  Comm ParallelismWithinCurrentPosses() const;
  virtual bool RemoveParallelization(Comm comm);
#endif //DOBLIS

  virtual const string& GetFunctionalityString() const;
  virtual bool IsShadow() const {return true;}
  virtual const RealPSet* GetReal() const {return m_realPSet;}
  virtual RealPSet* GetReal() {return m_realPSet;}
  virtual BasePSet* GetNewShadow();
};

