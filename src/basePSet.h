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
  bool m_currHasPrinted;
  bool m_hasProped;
  bool m_isTopLevel;
  BasePSet();
  virtual ~BasePSet() = 0;
  virtual GraphNum NumPosses() const = 0;
  virtual bool operator==(const BasePSet &rhs) const = 0;
  virtual void Prop() = 0;
  virtual void ClearBeforeProp();
  Node* InTun(unsigned int num) const;
  Node* OutTun(unsigned int num) const;
  bool GlobalSimplification(const TransMap &globalSimplifiers, const TransMap &simplifiers);
  virtual bool CanMerge(BasePSet *pset) const;
  virtual bool IsTransparent() const {return true;}
  virtual GraphNum TotalCount() const = 0;
  virtual void Duplicate(const BasePSet *orig, NodeMap &map, bool possMerging);
  virtual void ClearPrinted();
  virtual bool IsLoop() const {return false;}
  void RemoveInTun(Node *tun);
  void RemoveOutTun(Node *tun);
  void FormSetAround();
  bool CanPrint() const;
  virtual BasePSet* GetNewInst() = 0;

  void Flatten(ofstream &out) const;
  virtual void FlattenCore(ofstream &out) const {}
  void Unflatten(ifstream &in, SaveInfo &info);
  virtual void UnflattenCore(ifstream &in, SaveInfo &info) {}

  virtual void PrePrint(IndStream &out, Poss *poss) {}
  virtual void PostPrint(IndStream &out, Poss *poss) {}

  virtual void BuildDataTypeCache();
  virtual void ClearDataTypeCache();

  virtual const PossMMap& GetPosses() const = 0;
  virtual PossMMap& GetPosses() = 0;

#if DOBLIS
  virtual bool IsCritSect() const {return false;}
  Comm ParallelismWithinCurrentPosses() const;
  virtual bool RemoveParallelization(Comm comm) = 0;
#endif //DOBLIS

  virtual const string& GetFunctionalityString() const = 0;
  virtual bool IsReal() const {return false;}
  virtual bool IsShadow() const {return false;}
};


typedef vector<BasePSet*> PSetVec;
typedef PSetVec::iterator PSetVecIter;
typedef PSetVec::const_iterator PSetVecConstIter;
typedef set<const BasePSet*> PSetSet;
typedef PSetSet::iterator PSetSetIter;
typedef PSetSet::const_iterator PSetSetConstIter;
