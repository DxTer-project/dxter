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

#include "LLDLA.h"
#include "transform.h"

#if DOLLDLA

class SVMulPackResidualToVRW : public SingleTrans {
 private:
  Layer m_fromLayer, m_toLayer;
  VecType m_vecType;

 public:
 SVMulPackResidualToVRW(Layer fromLayer, Layer toLayer, VecType vecType)
   : m_fromLayer(fromLayer), m_toLayer(toLayer), m_vecType(vecType) {}

  virtual string GetType() const;
  virtual bool IsRef() const { return true; }

  virtual bool CanApply(const Node* node) const;
  virtual void Apply(Node* node) const;
};

#endif // DOLLDLA
