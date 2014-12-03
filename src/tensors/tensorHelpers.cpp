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


#include "layers.h"
#if DOTENSORS
#include "base.h"
#include "tensorHelpers.h"
#include "string.h"
#include "helperNodes.h"

using namespace std;

Copy::Copy(Layer layer)
{ 
  SetLayer(layer);
}

NodeType Copy::GetType() const
{
  return "Copy " + LayerNumToStr(GetLayer());
}

Phase Copy::MaxPhase() const 
{
  if (m_layer == DM1LAYER || m_layer == DM2LAYER)
    throw;
  else if (m_layer == SMLAYER)
    return NUMPHASES;
  else
    throw;
}

void Copy::PrintCode(IndStream &out)
{
  out.Indent();
  *out << GetInputName(0).str() << " = "
       << GetInputName(1).str() << ";\n";
}

void Copy::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    if (m_layer == ABSLAYER)
      m_cost = TotalNumberOfElements(0);
    else if (m_layer == DM1LAYER || m_layer == DM2LAYER)
      m_cost = TotalNumberOfElements(0);
    else if (m_layer == SMLAYER) {
      m_cost = TotalNumberOfLocalElements(0);
      if (InputDataType(0).GetEffectiveDist() != InputDataType(1).GetEffectiveDist())
	throw;

    }
    else {
      cout << LayerNumToStr(m_layer) << endl;
      throw;
    }
  }
}

#endif // DOTENSORS
