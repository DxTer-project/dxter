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
#include "intLoop.h"
#include "shadowPSet.h"

class ShadowLoop : public IntLoop<ShadowPSet>
{
  virtual BasePSet* GetNewInst() {return (BasePSet*)(new ShadowLoop);}
  virtual BasePSet* GetShadow();
  virtual int GetBS() const;
  virtual BSSize GetBSSize() const;
#if TWOD
  virtual DimName GetDimName() const;
#endif
  virtual void Prop();
  virtual LoopType GetType() const;
#if DOBLIS
  bool HasIndepIters() const;
  bool IsParallel() const {return m_comm!=CORECOMM;}
#endif
};
