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

#include <string>

#include "allTransformations.h"
#include "DLAReg.h"
#include "oneStageTimingResult.h"
#include "runtimeEvaluation.h"

#define TIMEANDCULLBEFOREUNROLLING 0

static string evalDirName = "runtimeEvaluation";
static SanityCheckSetting sanityCheckSetting = CHECKOUTPUTBUFFERS;
static TimingSetting timingSetting = TWOPHASETIMING;
static unsigned int numberOfImplementationsToEvaluate = 1000;
static int minCycles = 100000000;

#if TIMEANDCULLBEFOREUNROLLING
static double percentToKeep = .5;
void FirstPhaseTimingAndCulling(LLDLAUniverse* uni, ProblemInstance* problemInstance, double percentToCull);
#endif


ProblemInstanceStats* RunProblemWithRTE(int algNum, RealPSet* algPSet, ProblemInstance* problemInstance) {
  auto uni = RunProblem(algNum, algPSet, problemInstance);
  auto pStats = RuntimeEvaluation(algNum, uni, problemInstance);
  delete uni;
  return pStats;
}

std::string LLDLAPhaseString(LLDLAPhase phase) {
  switch(phase) {
  case(LLDLALOOPPHASE):
    return "LLDLA loop phase";
  case(LLDLALOOPUNROLLPHASE):
    return "LLDLA loop unroll phase";
  case(LLDLAPRIMPHASE):
    return "LLDLA primitive phase";
  case(LLDLARTLPHASE):
    return "LLDLA RTL phase";
  }
  LOG_FAIL("replacement for throw call");
  throw;
}

void RunPhase(Universe* uni, int numIters, LLDLAPhase phase) {
  time_t start, end;

  cout << "Expanding " << LLDLAPhaseString(phase) << endl;

  time(&start);
  uni->Expand(numIters, phase, LLDLACull);
  time(&end);

  cout << LLDLAPhaseString(phase) << " took " << difftime(end,start) << " seconds\n";
  cout << "Propagating\n";

  cout.flush();
  time(&start);
  uni->Prop();
  time(&end);

  cout << "Propagation took " << difftime(end,start) << " seconds\n";
}

LLDLAUniverse* RunProblem(int algNum, RealPSet* startSet, ProblemInstance* problemInstance) {
  RegAllLLDLANodes();
  AddTransformations();

  int numIters = -1;
  auto uni = new LLDLAUniverse();
  time_t start, end;

  uni->PrintStats();

  cout << "Initializing universe with start set" << endl;
  uni->Init(startSet);
  cout << "Initialized universe\n";

  cout << "Setting up problem" << endl;
  uni->SetUpOperation(startSet);
  cout << "Done with problem setup" << endl;

  problemInstance->SetCost(uni->GetOperationFlopCost());
  cout << "Implementation for correctness check:\n" << uni->GetSanityCheckImplStr();
  cout << "Flops for operation = " << std::to_string((long double) uni->GetOperationFlopCost()) << endl;

  time(&start);

  if ((CurrPhase == LLDLALOOPPHASE) && DOLLDLALOOPPHASE) {
    RunPhase(uni, numIters, LLDLALOOPPHASE);
  }

  if ((CurrPhase == LLDLARTLPHASE) && DOLLDLARTLPHASE) {
    RunPhase(uni, numIters, LLDLARTLPHASE);
  }

  if ((CurrPhase == LLDLAPRIMPHASE) && DOLLDLAPRIMPHASE) {
    RunPhase(uni, numIters, LLDLAPRIMPHASE);
  }

#if TIMEANDCULLBEFOREUNROLLING
  GraphNum num = uni->TotalCount();
  FirstPhaseTimingAndCulling(uni, problemInstance, percentToKeep);
  cout << "Aftering first timing phase, " << num << " -> " << uni->TotalCount() << " graphs left\n";
#endif

  if ((CurrPhase == LLDLALOOPUNROLLPHASE) && DOLLDLALOOPUNROLLPHASE) {
    RunPhase(uni, numIters, LLDLALOOPUNROLLPHASE);
  }

  time(&end);
  cout << "Full expansion took " << difftime(end,start) << " seconds\n";
  cout.flush();

  uni->ClearTransformations();

  return uni;
}

ProblemInstanceStats* RuntimeEvaluation(int algNum, LLDLAUniverse* uni, ProblemInstance* problemInstance) {
  LOG_A("Starting runtime evaluation for " + problemInstance->GetName());
  cout << "Writing all implementations to runtime eval files\n";
  RuntimeTest rtest(problemInstance, uni, minCycles);
  RuntimeEvaluator evaler = RuntimeEvaluator(evalDirName);

  cout << "About to evaluate\n";
  auto impMap = uni->ImpStrMap(false, numberOfImplementationsToEvaluate);
  vector<TimingResult*>* timingResults = evaler.EvaluateImplementations(sanityCheckSetting, timingSetting, rtest, impMap.get(), uni->GetSanityCheckImplStr());
  cout << "Done evaluating\n";

  vector<OneStageTimingResult*>* oneStageResults = reinterpret_cast<vector<OneStageTimingResult*>*>(timingResults);
  auto pStats = new ProblemInstanceStats(problemInstance, oneStageResults);
  pStats->PrettyPrintPerformanceStats();

  GraphNum best = pStats->GetBestAvgFlopsPerCycleImpl();
  cout << "Best Avg. flops/cycle = " << pStats->GetBestAvgFlopsPerCycle() << endl;

  cout << impMap.get()->find(best)->second.str << endl;

  LOG_A("Done with runtime evaluation of " + problemInstance->GetName());

  for (auto elem : *timingResults)
    delete elem;
  delete timingResults;

  return pStats;
}

#if TIMEANDCULLBEFOREUNROLLING
void FirstPhaseTimingAndCulling(LLDLAUniverse* uni, ProblemInstance* problemInstance, double percentToCull)
{
  LOG_A("Starting first-phase runtime evaluation for " + problemInstance->GetName());
  cout << "Writing all implementations to runtime eval files\n";
  RuntimeTest rtest(problemInstance, uni, minCycles);
  RuntimeEvaluator evaler = RuntimeEvaluator(evalDirName);

  cout << "About to evaluate\n";

  auto impMap = uni->ImpStrMap(true, numberOfImplementationsToEvaluate);
  vector<TimingResult*>* timingResults = evaler.EvaluateImplementations(sanityCheckSetting, timingSetting, rtest, impMap.get(), uni->GetSanityCheckImplStr());
  cout << "Done evaluating\n";

  vector<OneStageTimingResult*>* oneStageResults = reinterpret_cast<vector<OneStageTimingResult*>*>(timingResults);
  auto pStats = new ProblemInstanceStats(problemInstance, oneStageResults);

  vector<GraphNum> keepers;
  pStats->GetNBest(keepers,ceil(impMap->size()*percentToCull));

  uni->m_pset->ClearKeeperFromAll();

  cout << "should have " << ceil(impMap->size()*percentToCull) << endl;
  cout << "keepers size " << keepers.size() << endl;

  for(auto num : keepers) {
    (*impMap)[num].iter->SetCurrAsKeeper();
  }

  for (auto &info : *impMap) {
    delete info.second.iter;
  }


  for (auto elem : *timingResults)
    delete elem;
  delete timingResults;


  uni->m_pset->DeleteNonKeepers();


  
  LOG_A("Done with first-phase runtime evaluation of " + problemInstance->GetName());
}
#endif //TIMEANDCULLBEFOREUNROLLING

#endif // DOLLDLA
