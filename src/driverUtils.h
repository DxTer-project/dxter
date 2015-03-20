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

#ifndef DRIVER_UTILS_H_
#define DRIVER_UTILS_H_

#pragma once

#include "base.h"
#include "LLDLA.h"

Trans CharToTrans(char c);
Tri CharToTri(char c);
Side CharToSide(char c);

#if DOLLDLA
VecType CharToVecType(char c);

double BestFlopsPerCycle(Type type, ImplementationRuntimeMap &impTimes, double flopCost);
GraphNum PrintImpMapStats(Type type, ImplementationRuntimeMap &impTimes, double flopCost);

#endif // DOLLDLA

Type CharToType(char c);

#endif // DRIVER_UTILS_H_
