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

#include "mvmulPack.h"

#if DOLLDLA

#include "mvmul.h"
#include "pack.h"
#include "packingUtils.h"
#include "verticalUnpack.h"

bool MVMulPackOutput::CanApply(const Node* node) const {
  if (node->GetNodeClass() == MVMul::GetClass()) {
    const MVMul* mvmul = static_cast<const MVMul*>(node);
    if (mvmul->GetLayer() != m_fromLayer) {
      return false;
    }
    return !(mvmul->InputMIsMultipleOfVecRegWidth(0))
      && mvmul->GetInputM(0) > 0;
  }
  throw;
}

void MVMulPackOutput::Apply(Node* node) const {
  cout << endl << endl << "APPLYING MVMULPACKOUTPUT" << endl << endl;
  cout << endl << "Calling Prop on node" << endl;
  node->Prop();
  cout << "Done with node->Prop()" << endl << endl;
  //  cout << "Input node input 0 num rows = ";
  //  DLANode* dla = static_cast<DLANode*>(node);
  //  cout << dla->GetInputNumRows(0) << endl;
  Pack* packA = PackToMultipleOf(m_toLayer, node->Input(0), node->InputConnNum(0), node, 0, DIMM, node->GetVecRegWidth());

  Pack* packY = PackToMultipleOf(m_toLayer, node->Input(2), node->InputConnNum(2), node, 2, DIMM, node->GetVecRegWidth());

  MVMul* newMVMul = new MVMul(m_toLayer);
  newMVMul->AddInputs(6,
		      packA, 0,
		      node->Input(1), node->InputConnNum(1),
		      packY, 0);

  cout << "CHECKING THAT MVMULPACKOUTPUT DOESN'T INTRODUCE MIDSIZES" << endl;
  if (!newMVMul->InputDimsAreOneRepeatedSizeEach(0)
      || !newMVMul->InputDimsAreOneRepeatedSizeEach(1)
      || !newMVMul->InputDimsAreOneRepeatedSizeEach(2)) {
    cout << "ERROR IN THE MVMulPackOutput transformation" << endl;
    cout << "ERROR: MVMul input dimensions are not single, repeated sizes" << endl;
    throw;
  }
  cout << "DONE, MVMULPACKOUTPUT DID NOT INTRODUCE A MIDSIZE" << endl;

  Unpack* unpack = new VerticalUnpack(m_toLayer);
  unpack->AddInputs(4,
		    newMVMul, 0,
		    node->Input(2), node->InputConnNum(2));

  node->m_poss->AddUp(node->m_poss->m_possNodes, unpack, false, true);
  node->RedirectChildren(unpack, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#endif // DOLLDLA
