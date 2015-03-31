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

#ifndef RUN_BENCHMARK_H_
#define RUN_BENCHMARK_H_

#include "blasExamples.h"

#if DOLLDLA

void SVMulAddBenchmark(Type type, int m, int mInc, int numIters);
void SVMulBenchmark(Type type, VecType vecType, int m, int mInc, int iters);
void ColVAddBenchmark(Type type, int m, int increment, int numIters);
void MAddBenchmark(Type type, int mBase, int mInc, int nBase, int nInc, int numIters);
void DotProductBenchmark(Type type, int m, int increment, int numIters);
void RunDotProdBenchmarks();
void RunAxpyBenchmarks();
void RunGemvBenchmarks();
void RunAllBenchmarks();

#endif // DOLLDLA

#endif
