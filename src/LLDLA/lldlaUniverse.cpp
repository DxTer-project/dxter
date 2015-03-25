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

#include "lldlaUniverse.h"

#if DOLLDLA

#include <sstream>

#include "helperNodes.h"

void LLDLAUniverse::Init(RealPSet* seed) {
Universe::Init(seed);
  SetupFunctionArguments(seed);
}

void LLDLAUniverse::SetupFunctionArguments(RealPSet* seed) {
  int pSize = seed->m_posses.size();
  cout << std::to_string((long long int) pSize) << endl;
  Poss *poss = seed->m_posses.begin()->second;
  for (auto node : poss->m_possNodes) {
    if (node->GetNodeClass() == InputNode::GetClass()) {
      InputNode* inNode = (InputNode*) node;
      m_declarationVectors.push_back(inNode->DataDeclaration());
      m_constantDefines.push_back(inNode->NumRowsDefine());
      m_constantDefines.push_back(inNode->NumColsDefine());
      m_constantDefines.push_back(inNode->RowStrideDefine());
      m_constantDefines.push_back(inNode->ColStrideDefine());
      m_argNames.push_back(inNode->GetName(0).str());
    }
  }
}

void LLDLAUniverse::SetUpOperation(RealPSet* startSet) {
  this->Prop();
  GraphIter* graphIter = new GraphIter(startSet->m_posses.begin()->second);
  m_flopCost = graphIter->EvalAndSetBest();
  std::stringstream ss;
  IndStream optOut(&ss, LLDLASTREAM);
  graphIter->PrintRoot(optOut, 0, true, startSet);
  m_sanityCheckImplStr = ss.str();
}

#endif // DOLLDLA
