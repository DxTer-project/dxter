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

#include "copy.h"
#include "DLAReg.h"
#include "exampleUtils.h"
#include "mvmul.h"
#include "nodeTestExamples.h"
#include "pack.h"
#include "partition.h"
#include "recombine.h"
#include "setToZero.h"
#include "svmul.h"
#include "verticalUnpack.h"
#include "verticalPack.h"

#if DOLLDLA

RealPSet* VerticalPackUnpackTest(Type dataType, int m) {
  auto tunX = InputTunnel("x",
			  m - 3, 1,
			  1, m - 3,
			  dataType);

  auto tunY = InputTunnel("y",
			  m, 1,
			  1, m,
			  dataType);

  auto pack = new VerticalPack(ABSLAYER);
  pack->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  auto unpack = new VerticalUnpack(ABSLAYER);
  unpack->AddInputs(4,
		    pack, 0,
		    tunX, 0);

  return WrapInPSet(unpack);
}

RealPSet* VerticalRefinedPackTest(Type dataType, int m) {
  auto tunX = InputTunnel("x",
			  m - 3, 1,
			  1, m - 3,
			  dataType);

  auto tunY = InputTunnel("y",
			  m, 1,
			  1, m,
			  dataType);

  auto part = new Partition(ABSLAYER, VERTICAL, m - 3);
  part->AddInput(tunY, 0);

  auto copy = new Copy(ABSLAYER);
  copy->AddInputs(4,
		  tunX, 0,
		  part, 0);

  auto zero = new SetToZero(ABSLAYER);
  zero->AddInput(part, 1);

  auto recombine = new Recombine(ABSLAYER, VERTICAL);
  recombine->AddInputs(6,
		       copy, 0,
		       zero, 0,
		       tunY, 0);

  return WrapInPSet(recombine);
}

RealPSet* VerticalPartitionRecombineTest(Type dataType, int m) {
  auto tunX = InputTunnel("x",
			  1, m,
			  m, 1,
			  dataType);

  auto part = new Partition(ABSLAYER, HORIZONTAL, m - 5);
  part->AddInput(tunX, 0);

  auto rec = new Recombine(ABSLAYER, HORIZONTAL);
  rec->AddInputs(6,
		 part, 0,
		 part, 1,
		 tunX, 0);

  return WrapInPSet(rec);
}

RealPSet* HorizontalPartitionRecombineTest(Type dataType, int m) {
  auto xIn = new InputNode("x",
			   m, 1,
			   1, m,
			   dataType);

  auto tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  auto part = new Partition(ABSLAYER, VERTICAL, m - 5);
  part->AddInput(tunX, 0);

  auto rec = new Recombine(ABSLAYER, VERTICAL);
  rec->AddInputs(6,
		 part, 0,
		 part, 1,
		 tunX, 0);

  return WrapInPSet(rec);
}

RealPSet* CopyTest(Type dataType, int m, int n) {
  auto xIn = new InputNode("x",
			   m, n,
			   1, m,
			   dataType);

  auto yIn = new InputNode("y",
			   m, n,
			   1, m,
			   dataType);

  auto tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  auto tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  auto copy = new Copy(ABSLAYER);
  copy->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  Poss *innerPoss = new Poss(copy, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* PackTest(Type dataType, int m) {
  auto xIn = new InputNode("x",
			   m - 2, 1,
			   1, m - 2,
			   dataType);

  auto yIn = new InputNode("y",
			   m, 1,
			   1, m,
			   dataType);

  auto alpha = new InputNode("alpha",
			     1, 1,
			     1, 1,
			     dataType);

  auto tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn);

  auto tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn);

  auto tunAlpha = new Tunnel(POSSTUNIN);
  tunAlpha->AddInput(alpha);

  auto pack = new VerticalPack(ABSLAYER);
  pack->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  auto ax = new SVMul(COLVECTOR, ABSLAYER);
  ax->AddInputs(4,
		tunAlpha, 0,
		pack, 0);

  Poss *innerPoss = new Poss(ax, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* SetToZeroTest(Type dataType, int m, int n) {
  auto Ain = new InputNode("A",
			   n, m,
			   1, n,
			   dataType);

  auto xIn = new InputNode("x",
			   m, 1,
			   1, m,
			   dataType);

  auto yIn = new InputNode("y",
			   n, 1,
			   1, n,
			   dataType);

  auto tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  auto tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  auto tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  auto zeroY = new SetToZero(ABSLAYER);
  zeroY->AddInput(tunY, 0);

  MVMul* mvmul = new MVMul(ABSLAYER);
  mvmul->AddInputs(6,
		   tunA, 0,
		   tunX, 0,
		   zeroY, 0);

  return WrapInPSet(mvmul);
}

#endif // DOLLDLA
