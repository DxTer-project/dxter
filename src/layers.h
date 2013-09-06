/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2013, The University of Texas and Bryan Marker

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



#pragma once

#define DODM 1
#define DOSQM 0
#define DOSM 0

#if DODM

//Phases of generation
#define DPPHASE 0
#define ROPHASE 1

#define FIRSTPHASE DPPHASE

//Max phase for which to generate code
#define MAXPHASE ROPHASE
#define NUMPHASES (MAXPHASE+1)

#define DODPPHASE 1
#define DOROPHASE 1
#define DOSR1PHASE 0
#define DOSR2PHASE 0
#define DOSOPHASE 0
#define DOSMPPHASE 0

//Software layers
enum Layers {
  ABSLAYER = 0,
  DMLAYER,
  SMLAYER,
  BADLAYER,
  S1LAYER,
  S2LAYER,
  S3LAYER,
  SOLAYER,
};
#elif DOSQM

#define SR1PHASE 0
#define SR2PHASE 1
#define SR3PHASE 2
#define SOPHASE 3

#define FIRSTPHASE SR1PHASE

//Max phase for which to generate code
#define MAXPHASE SOPHASE
#define NUMPHASES (MAXPHASE+1)

#define DODPPHASE 0
#define DOROPHASE 0
#define DOSR1PHASE 1
#define DOSR2PHASE 1
#define DOSR3PHASE 1
#define DOSOPHASE 1
#define DOSMPPHASE 0

//Software layers
enum Layers {
  ABSLAYER = 0,
  S1LAYER,
  S2LAYER,
  S3LAYER,
  BADLAYER,
  DMLAYER,
  SMLAYER
};
#elif DOSM

#define SR1PHASE 0
#define SR2PHASE 1
#define SR3PHASE 2
#define SOPHASE 3
#define SMPPHASE 4

#define FIRSTPHASE SR1PHASE

//Max phase for which to generate code
#define MAXPHASE SMPPHASE
#define NUMPHASES (MAXPHASE+1)

#define DODPPHASE 0
#define DOROPHASE 0
#define DOSR1PHASE 1
#define DOSR2PHASE 1
#define DOSR3PHASE 1
#define DOSMPPHASE 1
#define DOSOPHASE 1

//Software layers
enum Layers {
  ABSLAYER = 0,
  S1LAYER,
  S2LAYER,
  S3LAYER,
  BADLAYER,
  DMLAYER,
  SMLAYER
};
#else
bad selection of layering
#endif

typedef unsigned int Layer;
