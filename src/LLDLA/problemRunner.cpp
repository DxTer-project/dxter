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

#include "allTransformations.h"
#include "DLAReg.h"
#include "runtimeEvaluation.h"

LLDLAUniverse* RunProblem(int algNum, RealPSet* algPSet, ProblemInstance* problemInstance) {
  RegAllLLDLANodes();
  AddTransformations();

  int numIters = -1;
  auto uni = new LLDLAUniverse();
  time_t start, start2, end;

  uni->PrintStats();

  cout << "Creating startSet\n";

  RealPSet *startSet = algPSet;
  
  cout << "Created startSet\n";

  uni->Init(startSet);
  
  cout << "Initialized universe\n";

  cout << "Setting up problem" << endl;

  uni->SetUpOperation(startSet);

  cout << "Done with problem setup" << endl;

  problemInstance->SetCost(uni->GetOperationFlopCost());
  cout << "IMPLEMENTATION FOR CORRECTNESS CHECK:\n" << uni->GetSanityCheckImplStr();
  cout << "Flops for operation = " << std::to_string((long double) uni->GetOperationFlopCost()) << endl;

  time(&start);

#if DOLLDLALOOPPHASE
  if (CurrPhase == LLDLALOOPPHASE) {

    cout << "Expanding LLDLA loop phase\n";

    uni->Expand(-1, LLDLALOOPPHASE, LLDLACull);
    time(&end);

    cout << "LLDLALOOP phase took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";

    cout.flush();
    time(&start2);
    uni->Prop();
    time(&end);

    cout << "Propagation took " << difftime(end,start2) << " seconds\n";

  }
#endif

#if DOLLDLALOOPUNROLLPHASE
  if (CurrPhase == LLDLALOOPUNROLLPHASE) {
    cout << "LLDLALOOPUNROLL phase\n";
    uni->Expand(-1, LLDLALOOPUNROLLPHASE, LLDLACull);
    time(&end);
    cout << "LLDLALOOPUNROLL phase took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni->Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOLLDLAPRIMPHASE
  if (CurrPhase == LLDLAPRIMPHASE) {
    cout << "Expanding LL DLA prim phase\n";
    cout << "Starting with " << uni->TotalCount() << endl;
    time(&start2);
    uni->Expand(numIters, LLDLAPRIMPHASE, LLDLACull);
    time(&end);
    cout << "LLDLAPRIM phase took " << difftime(end,start2) << " seconds\n";

    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni->Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

  cout << "Full expansion took " << difftime(end,start) << " seconds\n";
  cout.flush();

  return uni;
}

ProblemInstanceStats* RuntimeEvaluation(int algNum, LLDLAUniverse* uni, ProblemInstance* problemInstance) {
  cout << "Writing all implementations to runtime eval files\n";
  int minCycles = 100000000;
  RuntimeTest rtest(problemInstance->GetType(), problemInstance->GetName(), uni->m_argNames, uni->m_declarationVectors, uni->m_constantDefines, minCycles);
  string evalDirName = "runtimeEvaluation";
  RuntimeEvaluator evaler = RuntimeEvaluator(evalDirName);
  cout << "About to evaluate\n";
ImplementationRuntimeMap impMap = evaler.EvaluateImplementationsWithCorrectnessCheck(rtest, uni->ImpStrMap().get(), uni->GetSanityCheckImplStr());

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
