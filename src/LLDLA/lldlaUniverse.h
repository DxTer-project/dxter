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

#ifndef LLDLA_UNIVERSE_H_
#define LLDLA_UNIVERSE_H_

#include "universe.h"

#if DOLLDLA

class LLDLAUniverse : public Universe {
 private:
  string m_sanityCheckImplStr;
  Cost m_flopCost;

 public:
  vector<string> m_declarationVectors;
  vector<string> m_constantDefines;
  vector<string> m_argNames;

  LLDLAUniverse()
    : Universe::Universe() {}

  void Init(RealPSet* seed);
  void SetupFunctionArguments(RealPSet* seed);

  string GetSanityCheckImplStr() { return m_sanityCheckImplStr; }
  Cost GetOperationFlopCost() { return m_flopCost; }

  void SetUpOperation(RealPSet* startSet);
};

#endif // DOLLDLA

#endif
