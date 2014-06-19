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

#if DOLLDLA

#include <sstream>

#include "driverUtils.h"
#include "debug.h"
#include "LLDLAGemmTransformations.h"
#include "smmul.h"
#include "runtimeEvaluation.h"

Size one = 1;
Size smallSize = 10;
Size medSize = 100;
Size bigSize = 1000;
//Size bs = ELEM_BS;

PSet* GemmExample();
PSet* DoubleGemmExample();

Trans transA, transB;

ImplementationMap ImpStrMap(Universe *uni)
{
  ImplementationMap impMap;
  unsigned int i;
  for (i = 1; i <= uni->TotalCount(); i++) {
    std::stringbuf sbuf;
    std::ostream out(&sbuf);
    IndStream istream = IndStream(&out, LLDLASTREAM);
    uni->Print(istream, i);
    impMap.insert(NumImplementationPair(i, sbuf.str()));    
  }
  return impMap;
}

void AddTrans()
{
  //Introduces loops in the m, n, and k dimensions, respectively
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, USELLDLAMU), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, USELLDLAMU), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, USELLDLAMU), LLDLALOOPPHASE);

  //Introduces loops in the m and n dimension for SMMul
  Universe::AddTrans(SMMul::GetClass(), new SMulLoopRef(ABSLAYER, ABSLAYER, DIMM, USELLDLAMU), LLDLALOOPPHASE);
  Universe::AddTrans(SMMul::GetClass(), new SMulLoopRef(ABSLAYER, ABSLAYER, DIMN, USELLDLAMU), LLDLALOOPPHASE);

  //Lowers the layer tag of a Gemm node that is LLDLA_MU in all three dimensions
  Universe::AddTrans(Gemm::GetClass(), new LLDAGemmLowerLayer(ABSLAYER, LLDLAMIDLAYER, LLDLA_MU), LLDLALOOPPHASE);


  //Lowers the layer tag of a SMMul node that is LLDLA_MU in both dimensions
  Universe::AddTrans(SMMul::GetClass(), new SMulLowerLayer(ABSLAYER, LLDLAMIDLAYER, LLDLA_MU), LLDLALOOPPHASE);


  //Changes Gemm with transposition to non-transposed version use Transpose nodes
  Universe::AddTrans(Gemm::GetClass(), new GemmTransToNotTrans(LLDLAMIDLAYER), LLDLAPRIMPHASE);

  //Replaces a Gemm node with a LLDLAGemm node
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmToPrim(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER), LLDLAPRIMPHASE);

  //Lowers the layer tag of a SMMul node that is LLDLA_MU in both dimensions
  Universe::AddTrans(SMMul::GetClass(), new SMulLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER, LLDLA_MU), LLDLAPRIMPHASE);
}

void AddSimplifiers()
{ 
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
  PSet* (*algFunc)();
  //  unsigned int whichGraph = 0;
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
    uni.Init(algFunc());
    time(&start);
  }


#if DOLLDLALOOPPHASE
  if (CurrPhase == LLDLALOOPPHASE) {
    cout << "Expanding LL DLA loop phase\n";
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
  
  cout << "Writing all implementations to runtime eval files\n";
  
  string evalDirName = "runtimeEvaluation";
  string opName = "dxt_gemm";
  string prelude = "#include \"row_stride_lldla_primitives.h\"\n#define MUVALUE 2\nvoid dxt_gemm(int CNumCols, int CNumRows, int ANumCols,\ndouble *A, int ARowStride, double *B, int BRowStride, double *C, int CRowStride) {\ndouble b = 1.0;\ndouble *beta = &b;\ndouble *A1, *B1, *C1, *A11, *B11, *C11;\n";

  RuntimeEvaluator evaler = RuntimeEvaluator(evalDirName, opName, prelude);
  evaler.WriteImplementationsToFiles(ImpStrMap(&uni));
  
#if 1
  uni.PrintAll(algNum);
#else
  uni.PrintBest();
#endif

  /*  if (whichGraph <= 0)
      uni.PrintAll();
      else
      uni.Print(cout, CODE, whichGraph); */

  return 0;
}

PSet* GemmExample()
{
  InputNode *Ain = new InputNode("A input",  smallSize, smallSize, "A", 
				 NONUNITSTRIDE, UNITSTRIDE,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride");
  InputNode *Bin = new InputNode("B input",  smallSize, smallSize, "B", 
				 NONUNITSTRIDE, UNITSTRIDE,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride");
  InputNode *Cin = new InputNode("C input",  smallSize, smallSize, "C", 
				 NONUNITSTRIDE, UNITSTRIDE,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride");

  PossTunnel *tunA = new PossTunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  PossTunnel *tunB = new PossTunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  PossTunnel *tunC = new PossTunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  Gemm *gemm = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFBETA, REAL);
  gemm->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  tunC,0);

  Poss *innerPoss = new Poss(gemm,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
}

PSet* DoubleGemmExample()
{
  InputNode *Ain = new InputNode("A input",  smallSize, smallSize, "A", 
				 NONUNITSTRIDE, UNITSTRIDE,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride");
  InputNode *Bin = new InputNode("B input",  smallSize, smallSize, "B", 
				 NONUNITSTRIDE, UNITSTRIDE,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride");
  InputNode *Cin = new InputNode("C input",  smallSize, smallSize, "C", 
				 NONUNITSTRIDE, UNITSTRIDE,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride");

  PossTunnel *tunA = new PossTunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  PossTunnel *tunB = new PossTunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  PossTunnel *tunC = new PossTunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  Gemm *gemm1 = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFBETA, REAL);
  gemm1->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  tunC,0);

  Gemm *gemm2 = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFBETA, REAL);
  gemm2->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  gemm1,0);

  Poss *innerPoss = new Poss(gemm2,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
}





#endif //DOLLDLA
