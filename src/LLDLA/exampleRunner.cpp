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

#include "allTransformations.h"

#if DOLLDLA

#include <climits>
#include <sstream>
#include <time.h>

#include "exampleRunner.h"
#include "problemRunner.h"


ProblemInstanceStats* RunExample(int algNum, RealPSet* algPSet, ProblemInstance* problemInstance) {
  RegAllLLDLANodes();
  AddTransformations();

  int numIters = -1;
  Cost flopCost = 0;
  Universe uni;
  time_t start, start2, end;
  string absImpStr;

  uni.PrintStats();

  cout << "Creating startSet\n";

  RealPSet *startSet = algPSet;
  
  cout << "Created startSet\n";

  uni.Init(startSet);
  
  cout << "Initialized universe\n";
  
  uni.Prop();
  GraphIter* graphIter = new GraphIter(startSet->m_posses.begin()->second);
  cout << "Printing evaluation code\n";
  flopCost = graphIter->EvalAndSetBest();
  problemInstance->SetCost(flopCost);
  std::stringstream ss;
  IndStream optOut(&ss, LLDLASTREAM);
  graphIter->PrintRoot(optOut, 0, true, startSet);
  absImpStr = ss.str();

  cout << "IMPLEMENTATION FOR CORRECTNESS CHECK:\n" << absImpStr;
  cout << "Flops for operation = " << std::to_string((long double) flopCost) << endl;

  time(&start);

#if DOLLDLALOOPPHASE
  if (CurrPhase == LLDLALOOPPHASE) {

    cout << "Expanding LLDLA loop phase\n";

    uni.Expand(-1, LLDLALOOPPHASE, LLDLACull);
    time(&end);

    cout << "LLDLALOOP phase took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";

    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);

    cout << "Propagation took " << difftime(end,start2) << " seconds\n";

  }
#endif

#if DOLLDLALOOPUNROLLPHASE
  if (CurrPhase == LLDLALOOPUNROLLPHASE) {
    cout << "LLDLALOOPUNROLL phase\n";
    uni.Expand(-1, LLDLALOOPUNROLLPHASE, LLDLACull);
    time(&end);
    cout << "LLDLALOOPUNROLL phase took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOLLDLAPRIMPHASE
  if (CurrPhase == LLDLAPRIMPHASE) {
    cout << "Expanding LL DLA prim phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, LLDLAPRIMPHASE, LLDLACull);
    time(&end);
    cout << "LLDLAPRIM phase took " << difftime(end,start2) << " seconds\n";

    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

  cout << "Full expansion took " << difftime(end,start) << " seconds\n";
  cout.flush();

  auto pStats = RuntimeEvaluation(algNum, &uni, problemInstance, absImpStr);
  return pStats;
}

void RunExampleNoRuntimeEval(int algNum, RealPSet* algPSet, ProblemInstance* problemInstance) {
  RegAllLLDLANodes();
  AddTransformations();

  int numIters = -1;
  Cost flopCost = 0;
  Universe uni;
  time_t start, start2, end;
  string absImpStr;

  uni.PrintStats();

  cout << "Creating startSet\n";

  RealPSet *startSet = algPSet;
  
  cout << "Created startSet\n";

  uni.Init(startSet);
  
  cout << "Initialized universe\n";
  
  uni.Prop();
  GraphIter* graphIter = new GraphIter(startSet->m_posses.begin()->second);
  cout << "Printing evaluation code\n";
  flopCost = graphIter->EvalAndSetBest();
  problemInstance->SetCost(flopCost);
  std::stringstream ss;
  IndStream optOut(&ss, LLDLASTREAM);
  graphIter->PrintRoot(optOut, 0, true, startSet);
  absImpStr = ss.str();

  cout << "IMPLEMENTATION FOR CORRECTNESS CHECK:\n" << absImpStr;
  cout << "Flops for operation = " << std::to_string((long double) flopCost) << endl;

  time(&start);

#if DOLLDLALOOPPHASE
  if (CurrPhase == LLDLALOOPPHASE) {

    cout << "Expanding LLDLA loop phase\n";

    uni.Expand(-1, LLDLALOOPPHASE, LLDLACull);
    time(&end);

    cout << "LLDLALOOP phase took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";

    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);

    cout << "Propagation took " << difftime(end,start2) << " seconds\n";

  }
#endif

#if DOLLDLALOOPUNROLLPHASE
  if (CurrPhase == LLDLALOOPUNROLLPHASE) {
    cout << "LLDLALOOPUNROLL phase\n";
    uni.Expand(-1, LLDLALOOPUNROLLPHASE, LLDLACull);
    time(&end);
    cout << "LLDLALOOPUNROLL phase took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOLLDLAPRIMPHASE
  if (CurrPhase == LLDLAPRIMPHASE) {
    cout << "Expanding LL DLA prim phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, LLDLAPRIMPHASE, LLDLACull);
    time(&end);
    cout << "LLDLAPRIM phase took " << difftime(end,start2) << " seconds\n";

    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

  cout << "Full expansion took " << difftime(end,start) << " seconds\n";
  cout.flush();

  return;
}

#endif // DOLLDLA
