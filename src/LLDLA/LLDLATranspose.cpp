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

#include "LLDLATranspose.h"
#include "DLAOp.h"
#include "loopSupport.h"
#include "LLDLA.h"

#if DOLLDLA

LLDLATranspose::LLDLATranspose(Layer layer)
{
  m_layer = layer;
}

void LLDLATranspose::PrintCode(IndStream &out)
{
  DataTypeInfo info = InputDataType(0);
  out.Indent();
  *out << "int " << info.m_colStrideVar << "_tmp, " << info.m_numColsVar << "_tmp;\n";
  return;
}

Phase LLDLATranspose::MaxPhase() const
{
  switch (m_layer)
    { 
    case(ABSLAYER):
      return LLDLALOOPPHASE;
    case(LLDLAMIDLAYER):
      return LLDLAPRIMPHASE;
    case (LLDLAPRIMITIVELAYER):
      return NUMPHASES; 
    default:
      throw;
    }
}

NodeType LLDLATranspose::GetType() const
{
  return "LLDLATranspose";
}

Node* LLDLATranspose::BlankInst()
{
  return new LLDLATranspose(ABSLAYER);
}

void LLDLATranspose::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<1,1>::Duplicate(orig, shallow, possMerging);
  const LLDLATranspose *rhs = (LLDLATranspose*)orig;
  m_layer = rhs->m_layer;
  return;
}

void LLDLATranspose::Prop()
{

}

#endif // DOLLDLA
