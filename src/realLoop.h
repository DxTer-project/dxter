
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
#include "transform.h"
#include "basePSet.h"
#include "comm.h"
#include "intLoop.h"
#include "realPSet.h"

class RealLoop : public IntLoop<RealPSet>
{
 private:
#if TWOD
  DimName m_dim;
#endif

 public:
  //  typedef int Size;
  LoopType m_type;
  static int M_currLabel;
 public:
  IntSet m_label;
  BSSize m_bsSize;
#if DOBLIS
  Comm m_comm;
#endif
  int m_currIter;
  
  RealLoop();
  RealLoop(LoopType type);
  RealLoop(LoopType type, Poss *poss, BSSize bsSize);
  virtual ~RealLoop();
  virtual BasePSet* GetNewInst();
  virtual ShadowPSet* GetNewShadow();

  virtual const IntSet& GetLabel() const {return m_label;}
  virtual LoopType GetType() const {return m_type;}
  virtual void Duplicate(const BasePSet *orig, NodeMap &map, bool possMerging, bool useShadows);
  void AssignNewLabel();
  void SetBS(BSSize size);
  int GetBS() const;
  virtual BSSize GetBSSize() const {return m_bsSize;}
  void TryToDeleteLoopTunnelSetAndCleanUp(LoopTunnel *tun);
#if TWOD
  void SetDimName(DimName dim);
  DimName GetDimName() const {return m_dim;}
#endif

  void FillTunnelSizes();
  virtual void BuildDataTypeCache();

  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  static void FlattenStatic(ofstream &out);
  static void UnflattenStatic(ifstream &in);

#if DOBLIS
  void Parallelize(Comm comm);
  virtual bool HasIndepIters() const;
  bool IsParallel() const {return m_comm!=CORECOMM;}
#endif
  int GetCurrIter() const {return m_currIter;}
  void SetCurrIter(int iter) {m_currIter = iter;}
  inline bool IsUnrolled() const {return m_flags & SETLOOPISUNROLLED;}
  virtual Cost Prop();  
};
