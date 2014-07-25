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

enum BSSizeEnum { 
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
  USELLDLA2MU,
  USELLDLA3MU,
#endif
  USEUNITBS,
  BADBSSIZE 
};

class BSSize
{
 public:
  BSSizeEnum m_val;

 BSSize() : m_val(BADBSSIZE) {}

  explicit BSSize(BSSizeEnum val)
   : m_val(val) {}

  bool operator==(const BSSize &rhs) const {return m_val == rhs.m_val;}
  bool operator!=(const BSSize &rhs) const {return m_val != rhs.m_val;}

  string VarName() const;
  string SubSizeStr() const;
  string Str() const;

  inline Size Size() const {
    switch(m_val)
      {
#if DOELEM
      case (USEELEMBS):
	return ELEM_BS;
#elif DOBLIS
      case (USEBLISMC):
	return BLIS_MC_BS;
      case (USEBLISKC):
	return BLIS_KC_BS;
      case (USEBLISNC):
	return BLIS_NC_BS;
      case (USEBLISOUTERBS):
	return BLIS_OUTER_BS;
#elif DOTENSORS
      case (USETENSORBS):
	return TENSOR_BS;
#elif DOLLDLA
      case (USELLDLAMU):
	return LLDLA_MU;
      case (USELLDLA2MU):
	return 2*LLDLA_MU;
      case (USELLDLA3MU):
	return 3*LLDLA_MU;
#endif
      case (USEUNITBS):
	return ONE;
      default:
	throw;
      }
  }
};

#if DOELEM
extern BSSize ElemBS;
#elif DOBLIS
extern BSSize BlisMC;
extern BSSize BlisKC;
extern BSSize BlisNC;
extern BSSize BlisOuter;
#elif DOTENSORS
extern BSSize TensorBS;
#elif DOLLDLA
extern BSSize LLDLAMu;
extern BSSize LLDLA2Mu;
extern BSSize LLDLA3Mu;
#endif
extern BSSize BadBS;
extern BSSize UnitBS;


unsigned int GetNumElems(PartDir dir);
string PartDirToStr(PartDir dir);

class SplitBase;
class LoopTunnel;

class Loop : public PSet
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
  
  Loop();
  Loop(LoopType type);
  Loop(LoopType type, Poss *poss, BSSize bsSize);
  virtual PSet* GetNewInst() {return new Loop(m_type);}
  virtual bool IsLoop() const {return true;}
  virtual bool IsTransparent() const {return false;}
  virtual bool CanMerge(PSet *pset) const;
  virtual bool WorthFusing(Loop *loop);
  virtual void PrintCurrPoss(IndStream &out, GraphNum &graphNum);
  LoopType GetType() const {return m_type;}
  virtual void Duplicate(const PSet *orig, NodeMap &map, bool possMerging);
  void AssignNewLabel();
  void SetBS(BSSize size);
  int GetBS() const;
  SplitBase* GetControl() const;
  virtual void Prop();
  //  unsigned int NumIters() const;
  bool ValidIter() const;
  LoopTunnel* CreateNewLoopTunnels(Node *input, ConnNum num, Poss *possToCareAbout, UpStat stat);
  void TryToDeleteLoopTunnelSetAndCleanUp(LoopTunnel *tun);
#if TWOD
  void SetDimName(DimName dim);
  DimName GetDimName() const {return m_dim;}
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
