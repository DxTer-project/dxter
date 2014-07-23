/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2014, The University of Texas and Bryan Marker

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

#include "base.h"
#include "costs.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "transform.h"
#include "loopSupport.h"
#include <time.h>
#include "DLAReg.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "LLDLAGemm.h"
#include "vvdot.h"
#include <climits>


#if DOLLDLA

#define DOEMPIRICALEVAL 1
#define PRINTCOSTS 1

#define DOLOOPUNROLLING 1
#define DO2MUTRANSFORMATIONS 1
#define DO3MUTRANSFORMATIONS 1

#include <sstream>

#include "driverUtils.h"
#include "debug.h"
#include "LLDLAGemmTransformations.h"
#include "smmul.h"
#include "runtimeEvaluation.h"
#include "loopUnrolling.h"

Size one = 1;
Size smallSize = 12;
Size medSize = 36;
Size bigSize = 1000;
//Size bs = ELEM_BS;

RealPSet* GemmExample();
RealPSet* DoubleGemmExample();

Trans transA, transB;

ImplementationMap ImpStrMap(Universe *uni)
{
  ImplementationMap impMap;
  GraphNum i;
  for (i = 1; i <= uni->TotalCount(); i++) {
    std::stringbuf sbuf;
    std::ostream out(&sbuf);
    IndStream istream = IndStream(&out, LLDLASTREAM);
    uni->Print(istream, i);
    impMap.insert(NumImplementationPair(i, sbuf.str()));    
  }
  return impMap;
}

void PrintImpMap(ImplementationRuntimeMap &impTimes)
{
  ImplementationRuntimeMapIter mit;
  for (mit = impTimes.begin(); mit != impTimes.end(); ++mit) {
    TimeVecIter vit;
    cout << "IMPLEMENTATION # " << std::to_string(mit->first) << endl;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      cout << std::to_string(*vit) << endl;
    }
    cout << endl;
  }
}

GraphNum PrintImpMapInFlops(ImplementationRuntimeMap &impTimes, double flopCost, int chunkSize) {
  /***************************************************************************
   * WARNING: These numbers are processor specific to Dillon's machine in GDC
   ***************************************************************************/
  double ticksPerSec = 1.0e6;
  double peakFLOPS = 30e9;
  GraphNum bestImpNum = 0;
  double bestFLOPS = 0;
  ImplementationRuntimeMapIter mit;
  for (mit = impTimes.begin(); mit != impTimes.end(); ++mit) {
    TimeVecIter vit;
    cout << "IMPLEMENTATION # " << std::to_string(mit->first) << endl;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      double totalFlops = flopCost * chunkSize;
      double totalTimeInSecs = *vit / ticksPerSec;
      double actualFLOPS = totalFlops / totalTimeInSecs;
      double pctPeak = (actualFLOPS / peakFLOPS) * 100;
      if (actualFLOPS > bestFLOPS) {
	bestFLOPS = actualFLOPS;
	bestImpNum = mit->first;
      }
      cout << "FLOPS = " << std::to_string(actualFLOPS) << "\t%Peak = " << std::to_string(pctPeak) << endl;
      /*      if (pctPeak > 100) {
	cout << "pctPeak > 100\n";
	throw;
	}*/
    }
    cout << endl;
  }
  cout << "Best flops achieved: " << std::to_string(bestFLOPS) << endl;
  cout << "Best percent of peak: " << std::to_string((bestFLOPS / peakFLOPS) * 100) << endl;
  return bestImpNum;      
}

void AddTrans()
{
  //Introduces loops in the m, n, and k dimensions, respectively
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, USELLDLAMU), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, USELLDLAMU), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, USELLDLAMU), LLDLALOOPPHASE);

#if DO2MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, USELLDLA2MU), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, USELLDLA2MU), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, USELLDLA2MU), LLDLALOOPPHASE);
#endif 

#if DO3MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, USELLDLA3MU), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, USELLDLA3MU), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, USELLDLA3MU), LLDLALOOPPHASE);
#endif


  //Introduces loops in the m and n dimension for SMMul
  Universe::AddTrans(SMMul::GetClass(), new SMulLoopRef(ABSLAYER, ABSLAYER, DIMM, USELLDLAMU), LLDLALOOPPHASE);
  Universe::AddTrans(SMMul::GetClass(), new SMulLoopRef(ABSLAYER, ABSLAYER, DIMN, USELLDLAMU), LLDLALOOPPHASE);

  //Lowers the layer tag of a Gemm node that is USELLDLAMU in all three dimensions
  Universe::AddTrans(Gemm::GetClass(), new LLDAGemmLowerLayer(ABSLAYER, LLDLAMIDLAYER, BSSizeToSize(USELLDLAMU)), LLDLALOOPPHASE);


  //Lowers the layer tag of a SMMul node that is USELLDLAMU in both dimensions
  Universe::AddTrans(SMMul::GetClass(), new SMulLowerLayer(ABSLAYER, LLDLAMIDLAYER, USELLDLAMU), LLDLALOOPPHASE);

#if DOLOOPUNROLLING
#if DO2MUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(), 
  		     new FullyUnrollLoop(2), LLDLALOOPUNROLLPHASE);
#endif //DO2MUTRANSFORMATIONS

#if DO3MUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(), 
  		     new FullyUnrollLoop(3), LLDLALOOPUNROLLPHASE);
#endif //DO3MUTRANSFORMATIONS
#endif
  
  // Vector dot product transform
  Universe::AddTrans(VVDot::GetClass(), new VVDotLowerLayer(ABSLAYER, LLDLAMIDLAYER, USELLDLAMU), LLDLALOOPPHASE);
  
  Universe::AddTrans(VVDot::GetClass(), new VVDotLoopRef(ABSLAYER, ABSLAYER, USELLDLAMU), LLDLALOOPPHASE);
}

void AddSimplifiers()
{ 
  //Replaces a Gemm node with a LLDLAGemm node
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmToPrim(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER), SIMP);

  //Lowers the layer tag of a SMMul node that is USELLDLAMU in both dimensions
  Universe::AddTrans(SMMul::GetClass(), new SMulLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER, USELLDLAMU), SIMP);

  //Changes Gemm with transposition to non-transposed version use Transpose nodes
  Universe::AddTrans(Gemm::GetClass(), new GemmTransToNotTrans(LLDLAMIDLAYER), SIMP);
}

void Usage()
{
  cout << "./driver arg1 arg2 ...\n";
  cout <<" arg1 == 0  -> Load from file arg1\n";
  cout <<"         1  -> Gemm Example N/T N/T\n";
  cout <<"         2  -> Double Gemm Example N/T N/T\n";
}

int main(int argc, const char* argv[])
{
#ifdef _OPENMP
  omp_set_num_threads(1);
  omp_set_nested(true);
#endif

  //  PrintType printType = CODE;
  int numIters = -1;
  RealPSet* (*algFunc)();
  //  GraphNum whichGraph = 0;
  int algNum;
  string fileName;

  if(argc < 2) {
    Usage();
    return 0;
  }
  else {
    algNum = atoi(argv[1]);
    switch(algNum) {
    case(1):
      if (argc != 4) {
	Usage();
	return 0;
      }
      algFunc = GemmExample;
      transA = CharToTrans(*argv[2]);
      transB = CharToTrans(*argv[3]);
      break;
    case(2):
      if (argc != 4) {
	Usage();
	return 0;
      }
      algFunc = DoubleGemmExample;
      transA = CharToTrans(*argv[2]);
      transB = CharToTrans(*argv[3]);
      break;
    default:
      Usage();
      return 0;
    }
  }

  RegAllLLDLANodes();
  AddTrans();
  AddSimplifiers();

  Universe uni;
  time_t start, start2, end;
  uni.PrintStats();
  Cost flopCost = 0;
  string absImpStr;
  if (algNum==0) {
    time(&start);
    uni.Init(fileName);
    time(&end);
    cout << "Unflatten took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
  else {
    RealPSet *startSet = algFunc();
    uni.Init(startSet);
    uni.Prop();
    flopCost = startSet->EvalAndSetBest();
    // Print abstract implementation to string for use in testing
    // EXTREMELY HACKY, I could not figure out how to redirect an
    // ostream to a string
    std::stringstream ss;
    IndStream optOut(&ss, LLDLASTREAM);
    cout << "TEST\n";
    startSet->GetCurrPoss()->PrintRoot(optOut, 0, false);
    absImpStr = ss.str();
    cout << "IMPLEMENTATION FOR CORRECTNESS CHECK:\n" << absImpStr;
    cout << "Flops for operation = " << std::to_string(flopCost) << endl;
    time(&start);
  }


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

#if DOEMPIRICALEVAL  
  cout << "Writing all implementations to runtime eval files\n";

  int chunkSize = 3000;
  int numIterations = 1;
  RuntimeTest rtest("dxt_gemm", uni.m_argNames, uni.m_declarationVectors, uni.m_constantDefines, numIterations, chunkSize);
  string evalDirName = "runtimeEvaluation";
  RuntimeEvaluator evaler = RuntimeEvaluator(evalDirName);
  cout << "About to evaluate\n";
  ImplementationRuntimeMap impMap = evaler.EvaluateImplementationsWithCorrectnessCheck(rtest, ImpStrMap(&uni), absImpStr);
  cout << "Done evaluating\n";
  GraphNum best = PrintImpMapInFlops(impMap, flopCost, chunkSize);
  cout << "All implementations printed\n";
  cout << "Best times";

#endif //DOEMPIRICALEVAL

#if 1
  uni.PrintAll(algNum, best);
#else
  uni.PrintBest();
#endif

#if PRINTCOSTS  
  uni.PrintCosts(impMap);
#endif

  /*  if (whichGraph <= 0)
      uni.PrintAll();
      else
      uni.Print(cout, CODE, whichGraph); */

  return 0;
}

RealPSet* GemmExample()
{
  InputNode *Ain = new InputNode("A input", medSize, medSize, "A", 
				 medSize, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride");
  InputNode *Bin = new InputNode("B input", medSize, medSize, "B", 
				 medSize, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride");
  InputNode *Cin = new InputNode("C input",  medSize, medSize, "C", 
				 medSize, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride");

  PossTunnel *tunA = new PossTunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  PossTunnel *tunB = new PossTunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  PossTunnel *tunC = new PossTunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  Gemm *gemm = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, REAL);
  gemm->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  tunC,0);

  Poss *innerPoss = new Poss(gemm,true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* DoubleGemmExample()
{
  InputNode *Ain = new InputNode("A input",  medSize, smallSize, "A", 
				 smallSize, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride");
  InputNode *Bin = new InputNode("B input", smallSize, medSize, "B", 
				 medSize, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride");
  InputNode *Cin = new InputNode("C input",  medSize, medSize, "C", 
				 medSize, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride");

  PossTunnel *tunA = new PossTunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  PossTunnel *tunB = new PossTunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  PossTunnel *tunC = new PossTunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  Gemm *gemm1 = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, REAL);
  gemm1->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  tunC,0);

  Gemm *gemm2 = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, REAL);
  gemm2->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  gemm1,0);

  Poss *innerPoss = new Poss(gemm2,true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}


#endif //DOLLDLA

//  LocalWords:  DoubleGemmExample
