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

#include "packingUtils.h"

#if DOLLDLA

#include "horizontalPack.h"
#include "horizontalUnpack.h"
#include "localInput.h"
#include "madd.h"
#include "vadd.h"
#include "verticalPack.h"
#include "vvdot.h"
#include "verticalUnpack.h"

Node* CopySymmetricBinop(Node* binop) {
  Node* copy;
  if (binop->GetNodeClass() == MAdd::GetClass()) {
    copy = MAdd::BlankInst();
  } else if (binop->GetNodeClass() == VAdd::GetClass()) {
    copy = VAdd::BlankInst();
  } else if (binop->GetNodeClass() == VVDot::GetClass()) {
    copy = VVDot::BlankInst();
  } else {
    copy = VVDot::BlankInst();
    cout << "Error: " << binop->GetNodeClass() << " is not a symmetric, binary operation" << endl;
  }
  copy->Duplicate(binop, true, false);
  return copy;
}

int ComputePackedWidth(int length, int multiple) {
  if (length - (length % multiple) == 0) {
    return multiple;
  } else {
    return length - (length % multiple);
  }
}

Pack* PackToMultipleOf(Node* outNode, ConnNum outNum, Node* inNode, ConnNum inNum, DimName dim, int multiple) {
  /*  Pack* pack;
  LocalInput* locIn;
  DLANode* inN = (DLANode*) inNode;

  if (dim == DIMM) {
    int nRows = inN->GetInputNumRows(inNum);
    Size packDim = ComputePackedWidth(nRows, multiple);
    //    locIn = new LocalInput(inN->GetName(inNum) + std::to_string((unsigned long long) outNode),
			   
    pack = new VerticalPack(ABSLAYER);
  } else {
    int nCols = inN->GetInputNumCols(inNum);
    Size packDim = ComputePackedWidth(nCols, multiple);
    
    pack = new HorizontalPack(ABSLAYER);
  }

  pack->AddInputs(4,
		  outNode, outNum,
		  locIn, 0);
		  return pack;*/
  throw;
}

Unpack* PackBinarySymmetricOperation(Node* binop, DimName dim, int multiple) {
  auto operand0Pack = PackToMultipleOf(binop->Input(0), binop->InputConnNum(0), binop, 0, dim, multiple);
  auto operand1Pack = PackToMultipleOf(binop->Input(1), binop->InputConnNum(1), binop, 1, dim, multiple);

  auto mainBinop = CopySymmetricBinop(binop);
  mainBinop->AddInputs(4,
		       operand0Pack, 0,
		       operand1Pack, 0);

  auto residualBinop = CopySymmetricBinop(binop);
  residualBinop->AddInputs(4,
			   operand0Pack, 1,
			   operand1Pack, 1);

  Unpack* unpack;
  if (dim == DIMM) {
    unpack = new VerticalUnpack(ABSLAYER);
  } else {
    unpack = new HorizontalUnpack(ABSLAYER);
  }
  unpack->AddInputs(4,
		    mainBinop, 0,
		    residualBinop, 0);

  return unpack;
}

#endif // DOLLDLA
