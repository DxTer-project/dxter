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

#include "horizontalUnpack.h"
#include "madd.h"
#include "vadd.h"
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

Partition* PartitionIntoMainAndResidual(Node* node, ConnNum outNum, DimName dim, int multiple) {
  Partition* part;
  if (dim == DIMM) {
    part = new Partition(ABSLAYER, VERTICAL, 0);
    //    Size splitPoint = node->GetInput
  } else {
    part = new Partition(ABSLAYER, HORIZONTAL, 0);
  }
  part->AddInput(node, outNum);
  throw;
}

Unpack* PackBinarySymmetricOperation(Node* binop, DimName dim, int multiple) {
  auto operand0Pack = PartitionIntoMainAndResidual(binop->Input(0), binop->InputConnNum(0), dim, multiple);
  auto operand1Pack = PartitionIntoMainAndResidual(binop->Input(1), binop->InputConnNum(1), dim, multiple);

  auto newBinop = CopySymmetricBinop(binop);
  newBinop->AddInputs(4,
		      operand0Pack, 0,
		      operand1Pack, 0);

  Unpack* unpack;
  if (dim == DIMM) {
    unpack = new VerticalUnpack(ABSLAYER);
  } else {
    unpack = new HorizontalUnpack(ABSLAYER);
  }
  unpack->AddInput(newBinop, 0);

  return unpack;
}

#endif // DOLLDLA
