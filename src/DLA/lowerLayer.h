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


#include "layers.h"

class LowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
#if DOBLIS||DOLLDLA||DOELEM
  DimName m_dim;
  Size m_bs;
 LowerLayer(Layer fromLayer, Layer toLayer, DimName dim, Size bs)
   : m_fromLayer(fromLayer), m_toLayer(toLayer), m_dim(dim), m_bs(bs) {}
#elif DOTENSORS
  LowerLayer(Layer fromLayer, Layer toLayer)
    : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
#endif
  virtual bool IsRef() const {return true;}

};

