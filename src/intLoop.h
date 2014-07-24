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
#include "basePSet.h"
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
  USELLDLA2MU,
  USELLDLA3MU,
#endif
  USEUNITBS,
  BADBSSIZE 
};

Size BSSizeToSize(BSSize size);	      
string BSSizeToVarName(BSSize size);
string BSSizeToSubSizeStr(BSSize size);

unsigned int GetNumElems(PartDir dir);
string PartDirToStr(PartDir dir);

class SplitBase;
class LoopTunnel;

class LoopInterface
{
 public:
  virtual int GetBS() const = 0;
  virtual BSSize GetBSSize() const = 0;
  virtual LoopType GetType() const = 0;
#if TWOD
  virtual DimName GetDimName() const = 0;
#endif
  virtual bool IsLoop() const {return true;}
  virtual bool IsTransparent() const {return false;}
  virtual SplitBase* GetControl() const = 0;
  virtual bool CanMerge(BasePSet *pset) const = 0;
  virtual bool WorthFusing(BasePSet *pset) = 0;
  virtual unsigned int LoopLevel() const = 0;
#if DOBLIS
  virtual bool HasIndepIters() const = 0;
#endif
};

template <class PSetType>
class IntLoop : virtual public PSetType, public LoopInterface
{
 public:
  virtual bool CanMerge(BasePSet *pset) const;
  virtual bool WorthFusing(BasePSet *pset);
  SplitBase* GetControl() const;
  virtual void Prop();
  virtual unsigned int LoopLevel() const;
  virtual void PrePrint(IndStream &out, Poss *poss);
  virtual void PostPrint(IndStream &out, Poss *poss);
};
