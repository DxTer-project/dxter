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

#include "partition.h"

#if DOLLDLA

#include "pack.h"
#include "partition.h"
#include "recombine.h"
#include "unpack.h"

Partition* PartitionIntoMainAndResidual(Node* outNode, ConnNum outNum, Node* inNode, ConnNum inNum, DimName dim, int multiple);
Pack* PackToMultipleOf(Layer layer, Node* node, ConnNum outNum, DimName dim, int multiple);
Recombine* PartitionBinarySymmetricOperation(Node* node, ConnNum outNum, DimName dim, int multiple);
Unpack* PackBinarySymmetricOperation(Layer layer, Node* binop, DimName dim, int multiple);

#endif // DOLLDLA
