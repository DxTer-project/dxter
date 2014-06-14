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
#include "transform.h"
#include "pset.h"
#include "comm.h"

enum PartDir { PARTDOWN,
	       PARTRIGHT,
	       PARTDIAG,
	       PARTUPWARD,
	       PARTLEFT,
	       PARTDIAGBACK,
	       LASTPARTDIR };

enum MatPart { PART0, PART1, PART2,
	       PART00, PART01, PART02, PART10, PART11, PART12, PART20, PART21, PART22,
	       LASTPART };

enum LoopType { ELEMLOOP,
		BLISLOOP,
		LLDLALOOP,
		UNKNOWNLOOP };

enum UpStat { FULLUP,
	      PARTUP,
	      NOTUP,
	      BADUP };

enum Quad { TL, TR,
	     BL, BR,
	     LASTQUAD };

enum BSSize { 
#if DOELEM
  USEELEMBS,
#elif DOBLIS
  USEBLISMC,
  USEBLISKC,
  USEBLISNC,
  USEBLISOUTERBS,
#elif DOTENSORS
  USETENSORBS,
#elif DOLLDLA
  USELLDLAMU,
#endif
  BADBSSIZE 
};

Size BSSizeToSize(BSSize size);	      
string BSSizeToSubSizeStr(BSSize size);

unsigned int GetNumElems(PartDir dir);
string PartDirToStr(PartDir dir);

class Split;
class LoopTunnel;

class Loop : public PSet
{
  //  typedef int Size;
  LoopType m_type;
  static int M_currLabel;
 public:
  IntSet m_label;
  BSSize m_bsSize;
#if DOBLIS
  Comm m_comm;
#endif
#if TWOD
  DimName m_dim;
#endif
  
  Loop();
  Loop(LoopType type);
  Loop(LoopType type, Poss *poss, BSSize bsSize);
  virtual PSet* GetNewInst() {return new Loop(m_type);}
  virtual bool IsLoop() const {return true;}
  virtual bool IsTransparent() const {return false;}
  virtual bool CanMerge(PSet *pset) const;
  virtual bool WorthFusing(Loop *loop);
  virtual void PrintCurrPoss(IndStream &out, unsigned int &graphNum);
  LoopType GetType() const {return m_type;}
  virtual void Duplicate(const PSet *orig, NodeMap &map, bool possMerging);
  void AssignNewLabel();
  void SetBS(BSSize size);
  int GetBS() const;
  Split* GetControl() const;
  virtual void Prop();
  //  unsigned int NumIters() const;
  bool ValidIter() const;
  LoopTunnel* CreateNewLoopTunnels(Node *input, unsigned int num, Poss *possToCareAbout, UpStat stat);
  void TryToDeleteLoopTunnelSetAndCleanUp(LoopTunnel *tun);
#if TWOD
  inline void SetDimName(DimName dim) {m_dim = dim;}
#endif
  unsigned int LoopLevel() const;

  void FillTunnelSizes();
  virtual void BuildDataTypeCache();

  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  static void FlattenStatic(ofstream &out);
  static void UnflattenStatic(ifstream &in);

#if DOBLIS
  void Parallelize(Comm comm);
  bool HasIndepIters() const;
  bool IsParallel() const {return m_comm!=CORECOMM;}
  bool OnlyParallelizedOnNonIndependentData() const;
#endif
};
