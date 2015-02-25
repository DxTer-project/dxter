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

#include "base.h"
//#include "poss.h"
//#include "possTunnel.h"

#define USESHADOWS 0

class RealPSet;
class ShadowPSet;

#define SETHASPROPEDFLAG  (1L<<1)
#define SETTOPLEVELFLAG   (1L<<2)
#define SETLOOPISUNROLLED (1L<<3)
#define SETISDELETINGFLAG (1L<<3)
#define SETHASMIGRATED    (1L<<4)
#define SETCHECKEDFORDUP    (1L<<5)

unsigned int FindInTunVec(const TunVec &vec, const Tunnel *node);


class BasePSet
{
 public:
  TunVec m_inTuns;
  TunVec m_outTuns;
  Poss *m_ownerPoss;
  Flags m_flags;
  BasePSet();
  virtual ~BasePSet() {}
  virtual GraphNum NumPosses() const = 0;
  virtual bool operator==(const BasePSet &rhs) const = 0;
  virtual Cost Prop() = 0;
  virtual void ClearBeforeProp();
  Tunnel* InTun(unsigned int num) const;
  Tunnel* OutTun(unsigned int num) const;
  virtual bool CanMerge(BasePSet *pset) const;
  virtual bool IsTransparent() const {return true;}
  virtual GraphNum TotalCount() const = 0;
  virtual void Duplicate(const BasePSet *orig, NodeMap &map, bool possMerging, bool useShadows);
  virtual bool IsLoop() const {return false;}
  void RemoveInTun(Node *tun);
  void RemoveOutTun(Node *tun);
  void FormSetAround();
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

#if DOTENSORS
  void GetDistTypeSet(StrSet &set) const;
  bool CheckDistTypeSet(StrSet &set) const;
  bool HasRedist() const;
  bool HasPermutableOut(unsigned int tunNum) const;
  bool HasPermutableIn(unsigned int tunNum) const;
#endif
};



