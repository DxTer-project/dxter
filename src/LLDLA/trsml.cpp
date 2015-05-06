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

#include "trsml.h"

#if DOLLDLA

void TRSML::PrintCode(IndStream& out) {
  throw;
}

void TRSML::Prop() {
  throw;
}

Phase TRSML::MaxPhase() const {
  throw;
}

Node* TRSML::BlankInst() {
  throw;
}

void TRSML::Duplicate(const Node* orig, bool shallow, bool possMerging) {
  throw;
}

NodeType TRSML::GetType() const {
  return "TRSML";
}

#endif // DOLLDLA
