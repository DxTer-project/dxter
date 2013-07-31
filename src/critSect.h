/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2013, The University of Texas and Bryan Marker

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

#include "pset.h"

class CritSect : public PSet
{
 public:
 CritSect() : PSet() {}
 CritSect(Poss *poss) : PSet(poss) {}
  virtual PSet* GetNewInst() {return new CritSect;}
  virtual void PrintCurrPoss(IndStream &out, unsigned int &graphNum);
  virtual bool IsCritSect() const {return true;}
  virtual bool IsTransparent() const {return false;}
  virtual bool CanMerge(PSet *pset) const {return false;}
  virtual void SanityCheck();
  virtual void BuildSizeCache();
};

class CritSectTunnel : public PossTunnel
{
 public:
  Sizes *m_msizes, *m_nsizes;
  Sizes *m_mlsizes, *m_nlsizes;
  CritSectTunnel();
 CritSectTunnel(PossTunType type) : PossTunnel(type) {}
  ~CritSectTunnel();
  static Node* BlankInst() { return new CritSectTunnel;}
  virtual Node* GetNewInst() {return BlankInst(); }
  virtual NodeType GetType() const {return "CritSectTunnel";}
  virtual void SanityCheck();
  virtual unsigned int NumOutputs() const {return 1;}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "CritSectTunnel";}
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
  virtual void BuildSizeCache();
  virtual void ClearSizeCache();
  //  virtual PossTunnel* GetSetTunnel();
};
