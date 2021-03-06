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

#include <vector>
#include <memory>
#include "base.h"
#include "realPSet.h"
#include "linearization/graphIter.h"

struct ImplInfo {
  string str;
  GraphIter *iter;
};

typedef Node* (*ConstructorFunc)();
typedef std::map<GraphNum, ImplInfo> ImplementationMap;
typedef std::pair<GraphNum, ImplInfo> NumImplementationPair;
typedef map<ClassType,ConstructorFunc> ConsFuncMap;
typedef ConsFuncMap::iterator ConsFuncMapIter;
typedef ConsFuncMap::const_iterator ConsFuncMapConstIter;
typedef std::map<const Transformation*, int> SimpPhaseMap;
typedef map<ClassType,TransVec*> TransMap;
typedef TransMap::iterator TransMapIter;
typedef TransMap::const_iterator TransMapConstIter;


extern unsigned int CurrPhase;

//Holds on to all Transformations and generated graphs
//Does the generation, evaluation, printing, etc.
class Universe
{

 public:
  static TransMap M_trans[NUMPHASES];
  static TransMap M_simplifiers;
  static SimpPhaseMap M_simpPhaseMap;
  static TransPtrMap M_transNames;
  static TransNameMap M_transPtrs;
  bool TakeIter(unsigned int phase);
  RealPSet *m_pset;
  static unsigned int M_transCount[NUMPHASES+2];
  static ConsFuncMap M_consFuncMap;

  Universe();
  void Init(RealPSet *seed);
  void Init(string fileName);
  virtual ~Universe();
#if DOTENSORS
  void CheckMaxDims();
#endif

  GraphNum Expand(unsigned int numIters, unsigned int phase, CullFunction Cull);
  void EvalCosts(IndStream &out, GraphNum &graphNum);
  GraphIter EvalCostsAndSetBest(Cost &best);
  void Print(IndStream &out, GraphNum &graphNum, bool currOnly = false);
  static void AddTrans(const ClassType &classType, Transformation *trans, int phase);
  static void AddSimp(const ClassType &classType, Transformation *trans, int phase);
  static void AddToMaps(Transformation *trans);
  void Cull();
  void Prop();
  int PrintAll(int algNum, GraphNum optGraph = 0);
  void PrintBest();
  void PrintCosts(const ImplementationRuntimeMap &impTimes);
  GraphNum TotalCount() const;
  void ClearFullyExpanded();
  void PrintStats();
  void Simplify();

  void ClearTransformations();

  unique_ptr<ImplementationMap> ImpStrMap(bool includeIters, unsigned int numGraphs = 0);

  static void RegCons(ClassType type, ConstructorFunc func);
  static Node* GetBlankClassInst(ClassType type);

  void SaveToFile(string fileName) const;
  void Flatten(ofstream &out) const;
  void LoadFromFile(string fileName);
  void Unflatten(ifstream &in);
  void CullWorstPerformers(double percentToCull, int ignoreThreshold);
  void CullAllBut(int num);
  void InlineAllSets();
  void EnforceMemConstraint(Cost maxMem);
};
