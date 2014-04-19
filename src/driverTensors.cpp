/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2013, The University of Texas and Bryan Marker

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
#include "omp.h"
#include "transform.h"
#include "loopSupport.h"
#include <time.h>
#include "DLAReg.h"
#include <omp.h>
#include "contraction.h"
#include "tensorRedist.h"

#if DOTENSORS

#include "debug.h"

Size one = 1;
Size smallSize = 10;
Size medSize = 100;
Size bigSize = 1000;
//Size bs = ELEM_BS;

PSet* Cont1Example();
PSet* MartinsExample();

void AddTrans()
{
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatC(DMLAYER, SMLAYER), DPTENSORPHASE);
  
#if 0
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatAAllReduce(DMLAYER, SMLAYER), DPTENSORPHASE);
#endif
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatASumScatter(DMLAYER, SMLAYER), DPTENSORPHASE);
  //    Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatBAllReduce(DMLAYER, SMLAYER), DPTENSORPHASE);
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatBSumScatter(DMLAYER, SMLAYER), DPTENSORPHASE);
  
  Universe::AddTrans(SumScatterUpdateNode::GetClass(), new SeparateRedistFromSumScatter, SUMSCATTERTENSORPHASE);
  Universe::AddTrans(SumScatterUpdateNode::GetClass(), new MoveSumScatterRedistAfter, SUMSCATTERTENSORPHASE);
  
#if 1
  for(Dim dim = 0; dim < NUM_GRID_DIMS; ++dim) {
    Universe::AddTrans(RedistNode::GetClass(), new SplitRedistribs(dim), ROTENSORPHASE);
    Universe::AddTrans(RedistNode::GetClass(), new SingleIndexAllToAll(dim), ROTENSORPHASE);
    Universe::AddTrans(SumScatterUpdateNode::GetClass(), new SplitSumScatter(dim), SUMSCATTERTENSORPHASE);
  }
#endif
}

void AddSimplifiers()
{ 
   Universe::AddTrans(RedistNode::GetClass(), new RemoveNOPRedistribs, SIMP);
   Universe::AddTrans(RedistNode::GetClass(), new RemoveWastedRedist, SIMP);
   Universe::AddTrans(RedistNode::GetClass(), new CombineRedistribs, SIMP);
}

void Usage()
{
  cout << "./driver arg1 arg2 arg3 arg4\n";
  cout <<" arg1 == 0  -> Load from file arg1\n";
  cout <<"         1  -> Contraction (abcd,cdef,abef)\n";
  cout <<"         2  -> Martin's Example\n";
}

int main(int argc, const char* argv[])
{
  omp_set_num_threads(1);
  omp_set_nested(true);
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
      algFunc = Cont1Example;
      break;
    case(2):
      algFunc = MartinsExample;
      break;
    default:
      Usage();
      return 0;
    }
  }

  RegAllTensorNodes();
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
    //    uni.SanityCheck();
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


#if DODPTENSORPHASE
  if (CurrPhase == DPTENSORPHASE) {
    cout << "Expanding DP phase\n";
    uni.Expand(-1, DPTENSORPHASE, TenCullDP);
    time(&end);
    cout << "DP phase took " << difftime(end,start) << " seconds\n";

    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOROTENSORPHASE
  if (CurrPhase == SUMSCATTERTENSORPHASE) {
    cout << "SumScatterOpt phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, SUMSCATTERTENSORPHASE, TenCullRO);
    time(&end);
    cout << "SumScatter phase took " << difftime(end,start2) << " seconds\n";
    
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOROTENSORPHASE
  if (CurrPhase == ROTENSORPHASE) {
    cout << "Expanding RO phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, ROTENSORPHASE, TenCullRO);
    time(&end);
    cout << "RO phase took " << difftime(end,start2) << " seconds\n";
    
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

#if 0
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

PSet* Cont1Example()
{
  Sizes sizes[4];

  for (Dim dim = 0; dim < 4; ++dim)
    sizes[dim].AddRepeatedSizes(bigSize, 1, 1);

  InputNode *Ain = new InputNode("A input",  sizes, "A", "acd");
  InputNode *Bin = new InputNode("B input",  sizes, "B", "cefd");
  InputNode *Cin = new InputNode("C input",  sizes, "C", "aef");

  PossTunnel *tunA = new PossTunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  PossTunnel *tunB = new PossTunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  PossTunnel *tunC = new PossTunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  Contraction *cont = new Contraction(DMLAYER,COEFONE,COEFONE,REAL,(string)"cd");
  cont->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  tunC,0);

  Poss *innerPoss = new Poss(cont,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;

}

PSet* MartinsExample()
{
  Sizes sizes[4];

  //a-d = medium
  //i-l = big

  for (Dim dim = 0; dim < 4; ++dim)
    sizes[dim].AddRepeatedSizes(medSize, 1, 1);

  InputNode *Uin = new InputNode("U input",  sizes, "U", "abcd");

  sizes[2].ClearSizes();
  sizes[2].AddRepeatedSizes(bigSize,1,1);
  sizes[3].ClearSizes();
  sizes[3].AddRepeatedSizes(bigSize,1,1);
  
  InputNode *Vin = new InputNode("V input",  sizes, "V", "acik");
  InputNode *T1in = new InputNode("T1 input",  sizes, "T1", "cdij");
  InputNode *T2in = new InputNode("T2 input",  sizes, "T2", "bcjk");
  InputNode *T3in = new InputNode("T3 input",  sizes, "T3", "abkl");
  InputNode *T4in = new InputNode("T4 input",  sizes, "T4", "abij");


  sizes[0].ClearSizes();
  sizes[0].AddRepeatedSizes(medSize,1,1);
  sizes[1].ClearSizes();
  sizes[1].AddRepeatedSizes(medSize,1,1);

  InputNode *Win = new InputNode("W input",  sizes, "W", "ijkl");
  

  Sizes ones[2];

  for (Dim dim = 0; dim < 2; ++dim)
    ones[dim].AddRepeatedSizes(one, 1, 1);

  DistType epDist;
  epDist.SetToScalarNoRep();

  InputNode *epIn = new InputNode("ep input",  ones, epDist, "epsilon", "");
  //InputNode *epIn = new InputNode("ep input",  ones, "epsilon", "");

  InputNode *tempIn = new InputNode("Temp input",  sizes, "Accum", "abij");

  PossTunnel *tunU = new PossTunnel(POSSTUNIN);
  tunU->AddInput(Uin,0);

  PossTunnel *tunV = new PossTunnel(POSSTUNIN);
  tunV->AddInput(Vin,0);

  PossTunnel *tunW = new PossTunnel(POSSTUNIN);
  tunW->AddInput(Win,0);

  PossTunnel *tunT1 = new PossTunnel(POSSTUNIN);
  tunT1->AddInput(T1in,0);

  PossTunnel *tunT2 = new PossTunnel(POSSTUNIN);
  tunT2->AddInput(T2in,0);

  PossTunnel *tunT3 = new PossTunnel(POSSTUNIN);
  tunT3->AddInput(T3in,0);

  PossTunnel *tunT4 = new PossTunnel(POSSTUNIN);
  tunT4->AddInput(T4in,0);

  PossTunnel *tunIn = new PossTunnel(POSSTUNIN);
  tunIn->AddInput(tempIn,0);

  PossTunnel *tunOutVal = new PossTunnel(POSSTUNIN);
  tunOutVal->AddInput(epIn,0);

  Contraction *cont1 = new Contraction(DMLAYER,COEFONE,COEFZERO,REAL,(string)"cd");
  cont1->AddInputs(6,
		  tunU,0,
		  tunT1,0,
		  tunIn,0);


  Contraction *cont2 = new Contraction(DMLAYER,COEFONE,COEFONE,REAL,(string)"ck");
  cont2->AddInputs(6,
		   tunV,0,
		   tunT2,0,
		   cont1,0);


  Contraction *cont3 = new Contraction(DMLAYER,COEFONE,COEFONE,REAL,(string)"kl");
  cont3->AddInputs(6,
		   tunW,0,
		   tunT3,0,
		   cont2,0);

  Contraction *cont4 = new Contraction(DMLAYER,COEFONE,COEFZERO,REAL,(string)"abij");
  cont4->AddInputs(6,
		   tunT4,0,
		   cont3,0,
		   tunOutVal,0);

  Poss *innerPoss = new Poss(cont4,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *out = new OutputNode("output");
  out->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(out,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;

}

#endif //DOTENSORS
