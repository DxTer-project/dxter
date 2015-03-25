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

#include "problemRunner.h"

#if DOLLDLA

#include "runtimeEvaluation.h"

ProblemInstanceStats* RuntimeEvaluation(int algNum, Universe* uni, ProblemInstance* problemInstance, string absImpStr) {
  cout << "Writing all implementations to runtime eval files\n";
  int minCycles = 100000000;
  RuntimeTest rtest(problemInstance->GetType(), problemInstance->GetName(), uni->m_argNames, uni->m_declarationVectors, uni->m_constantDefines, minCycles);
  string evalDirName = "runtimeEvaluation";
  RuntimeEvaluator evaler = RuntimeEvaluator(evalDirName);
  cout << "About to evaluate\n";
  ImplementationRuntimeMap impMap = evaler.EvaluateImplementationsWithCorrectnessCheck(rtest, uni->ImpStrMap().get(), absImpStr);

  cout << "Done evaluating\n";
  ProblemInstanceStats* pStats = new ProblemInstanceStats(problemInstance, &impMap);
  pStats->PrettyPrintPerformanceStats();

  GraphNum best = pStats->GetBestAvgFlopsPerCycleImpl();
  cout << "Best Avg. flops/cycle = " << pStats->GetBestAvgFlopsPerCycle() << endl;

#if 1
  uni->PrintAll(algNum, best);
#else
  uni->PrintBest();
#endif

#if PRINTCOSTS
  uni->PrintCosts(impMap);
#endif

  return pStats;
}

#endif // DOLLDLA
