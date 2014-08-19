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
//#include "poss.h"
//#include "possTunnel.h"

#define USESHADOWS DOTENSORS

class RealPSet;
class ShadowPSet;

#define SETHASPROPEDFLAG     (1L<<1)
#define SETTOPLEVELFLAG   (1L<<2)

class BasePSet
{
 public:
  NodeVec m_inTuns;
  NodeVec m_outTuns;
  Poss *m_ownerPoss;
  Flags m_flags;
  BasePSet();
  virtual ~BasePSet() {}
  virtual GraphNum NumPosses() const = 0;
  virtual bool operator==(const BasePSet &rhs) const = 0;
  virtual Cost Prop() = 0;
  virtual void ClearBeforeProp();
  Node* InTun(unsigned int num) const;
  Node* OutTun(unsigned int num) const;
  virtual bool CanMerge(BasePSet *pset) const;
  virtual bool IsTransparent() const {return true;}
  virtual GraphNum TotalCount() const = 0;
  virtual void Duplicate(const BasePSet *orig, NodeMap &map, bool possMerging, bool useShadows);
  virtual bool IsLoop() const {return false;}
  void RemoveInTun(Node *tun);
  void RemoveOutTun(Node *tun);
  void FormSetAround();
  bool CanPrint(const GraphIter *graphIter) const;
  virtual BasePSet* GetNewInst() = 0;
  virtual ShadowPSet* GetNewShadow() = 0;
  inline bool IsTopLevel() const { return m_flags & SETTOPLEVELFLAG; }

  void Flatten(ofstream &out) const;
  virtual void FlattenCore(ofstream &out) const {}
  void Unflatten(ifstream &in, SaveInfo &info);
  virtual void UnflattenCore(ifstream &in, SaveInfo &info) {}

  virtual void PrePrint(IndStream &out, Poss *poss) {}
  virtual void PostPrint(IndStream &out, Poss *poss) {}

  virtual void BuildDataTypeCache() = 0;
  virtual void ClearDataTypeCache() = 0;

  virtual const PossMMap& GetPosses() const = 0;
  virtual PossMMap& GetPosses() = 0;

#if DOBLIS
  virtual bool IsCritSect() const {return false;}
  Comm ParallelismWithinCurrentPosses() const;
  virtual bool RemoveParallelization(Comm comm) = 0;
#endif //DOBLIS

  virtual void PatchAfterDuplicate(NodeMap &map) = 0;
  virtual const string& GetFunctionalityString() const = 0;
  virtual bool IsReal() const {return false;}
  virtual bool IsShadow() const {return false;}
  virtual const RealPSet* GetReal() const = 0;
  virtual RealPSet* GetReal() = 0;
};



