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



#pragma once

#include "base.h"
#include "poss.h"
#include "basePSet.h"
#include "linearize.h"

class GraphIter;
typedef GraphIter* GraphIterPtr;
//typedef vector<GraphIter*> GraphIterVec;
//typedef GraphIterVec::iterator GraphIterVecIter;


class GraphIter
{
 public:
  Poss *m_poss;
  PossMMapIter *m_setIters;
  GraphIterPtr *m_subIters;
  Cost m_cost;

  GraphIter(Poss *poss);
  GraphIter(const GraphIter &iter);
  ~GraphIter();
  void Init(Poss *poss);
  bool Increment();
  GraphIter& operator=(const GraphIter &rhs);

  void GetCurrTransVec(TransVec &transVec);
  void AddCurrPossVars(VarSet &set) const;

  void EvalRoot(IndStream &out, GraphNum &graphNum, GraphNum whichGraph, GraphNum &optGraph, Cost &optCost);
  Cost Eval(TransConstVec &transList);
  Cost Eval();
  Cost EvalAndSetBest();

  void PrintRoot(IndStream &out, GraphNum whichGraph, bool currOnly, BasePSet *owner);
  void Print(IndStream &out, BasePSet *owner, StrSet liveSet, Cost currCost);
};
