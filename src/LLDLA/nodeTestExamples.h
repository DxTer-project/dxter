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

#include "LLDLA.h"

#if DOLLDLA

RealPSet* TwoDHorizontalUnpackTest(Type dataType, int m, int n);
RealPSet* TwoDVerticalUnpackTest(Type dataType, int m, int n);
RealPSet* TwoDVerticalPackUnpackTest(Type dataType, int m, int n);
RealPSet* VerticalPackUnpackTest(Type dataType, int m);
RealPSet* VerticalRefinedPackTest(Type dataType, int m);
RealPSet* CopyTest(Type dataType, int m, int n);
RealPSet* HorizontalCopyTest(Type dataType, int m, int n);
RealPSet* VerticalPartitionRecombineTest(Type dataType, int m);
RealPSet* HorizontalPartitionRecombineTest(Type dataType, int m);
RealPSet* PackTest(Type dataType, int m);
RealPSet* SetToZeroTest(Type dataType, int m, int n);

#endif // DOLLDLA
