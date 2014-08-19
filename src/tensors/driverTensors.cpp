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
#include "contraction.h"
#include "tensorRedist.h"
#include "tensorPermute.h"

#if DOTENSORS

#include "debug.h"

Size one = 1;
Size smallSize = 10;
Size medSize = 100;
Size bigSize = 1000;
//Size bs = ELEM_BS;

RealPSet* Cont1Example();
RealPSet* MartinsExample();
RealPSet* MartinsExample2();

void AddTrans()
{
  MultiTrans *trans = new MultiTrans;
  trans->AddTrans(new DistContToLocalContStatC(DMLAYER, SMLAYER));
  trans->AddTrans(new DistContToLocalContStatASumScatter(DMLAYER, SMLAYER));
  trans->AddTrans(new DistContToLocalContStatBSumScatter(DMLAYER, SMLAYER));
  Universe::AddTrans(Contraction::GetClass(), trans, DPTENSORPHASE);
  /*
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatC(DMLAYER, SMLAYER), DPTENSORPHASE);
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatASumScatter(DMLAYER, SMLAYER), DPTENSORPHASE);
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatBSumScatter(DMLAYER, SMLAYER), DPTENSORPHASE);  
  */
#if 0
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatAAllReduce(DMLAYER, SMLAYER), DPTENSORPHASE);
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatBAllReduce(DMLAYER, SMLAYER), DPTENSORPHASE);
#endif

  
  Universe::AddTrans(SumScatterUpdateNode::GetClass(), new SeparateRedistFromSumScatter, SUMSCATTERTENSORPHASE);
  Universe::AddTrans(SumScatterUpdateNode::GetClass(), new MoveSumScatterRedistAfter, SUMSCATTERTENSORPHASE);
  
#if 1
  for(Dim dim = 0; dim < NUM_GRID_DIMS; ++dim) {
    Universe::AddTrans(RedistNode::GetClass(), new SplitRedistribs(dim), ROTENSORPHASE);
    Universe::AddTrans(RedistNode::GetClass(), new SingleIndexAllToAll(dim), ROTENSORPHASE);
    Universe::AddTrans(RedistNode::GetClass(), new DoubleIndexAllToAll(dim), ROTENSORPHASE);
    Universe::AddTrans(RedistNode::GetClass(), new SplitAllGathers(dim), ROTENSORPHASE);
    Universe::AddTrans(SumScatterUpdateNode::GetClass(), new SplitSumScatter(dim), SUMSCATTERTENSORPHASE);

    for(Dim dim2 = 0; dim2 < NUM_GRID_DIMS; ++dim2) {
      if (dim2 != dim) {
	//	Universe::AddTrans(RedistNode::GetClass(), new CombineDisappearingModes(dim, dim2), ROTENSORPHASE);
	Universe::AddTrans(RedistNode::GetClass(), new PermuteDistribution(dim, dim2), ROTENSORPHASE);
      }
    }
  }
#endif
}

void AddSimplifiers()
{ 
   Universe::AddTrans(RedistNode::GetClass(), new RemoveNOPRedistribs, SIMP);
   Universe::AddTrans(RedistNode::GetClass(), new RemoveWastedRedist, SIMP);
   Universe::AddTrans(RedistNode::GetClass(), new CombineRedistribs, SIMP);
   Universe::AddTrans(Permute::GetClass(), new LowerPermute, SIMP);
}

void Usage()
{
  cout << "./driver arg1 arg2 arg3 arg4\n";
  cout <<" arg1 == 0  -> Load from file arg1\n";
  cout <<"         1  -> Contraction (abcd,cdef,abef)\n";
  cout <<"         2  -> Martin's Example\n";
  cout <<"         3  -> Martin's Example - real\n";
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
      algFunc = Cont1Example;
      break;
    case(2):
      algFunc = MartinsExample;
      break;
    case(3):
      algFunc = MartinsExample2;
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

#if DOSUMSCATTERTENSORPHASE
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

RealPSet* Cont1Example()
{
  Sizes sizes[4];

  for (Dim dim = 0; dim < 4; ++dim)
    sizes[dim].AddRepeatedSizes(bigSize, 1, 1);

  InputNode *Ain = new InputNode("A input",  sizes, "A", 3);
  InputNode *Bin = new InputNode("B input",  sizes, "B", 4);
  InputNode *Cin = new InputNode("C input",  sizes, "C", 3);

  Tunnel *tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  Tunnel *tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  Tunnel *tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  Contraction *cont = new Contraction(DMLAYER,COEFONE,COEFONE,REAL,"acd", "cefd", "aef", (string)"cd");
  cont->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  tunC,0);

  Poss *innerPoss = new Poss(cont,true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;

}

RealPSet* MartinsExample()
{
  Sizes sizes[4];

  //a-d = medium
  //i-l = big

  for (Dim dim = 0; dim < 4; ++dim)
    sizes[dim].AddRepeatedSizes(medSize, 1, 1);

  InputNode *Uin = new InputNode("U input",  sizes, "U", 4);

  sizes[2].ClearSizes();
  sizes[2].AddRepeatedSizes(bigSize,1,1);
  sizes[3].ClearSizes();
  sizes[3].AddRepeatedSizes(bigSize,1,1);
  
  InputNode *Vin = new InputNode("V input",  sizes, "V", 4);
  InputNode *T1in = new InputNode("T1 input",  sizes, "T1", 4);
  InputNode *T2in = new InputNode("T2 input",  sizes, "T2", 4);
  InputNode *T3in = new InputNode("T3 input",  sizes, "T3", 4);
  InputNode *T4in = new InputNode("T4 input",  sizes, "T4", 4);


  sizes[0].ClearSizes();
  sizes[0].AddRepeatedSizes(medSize,1,1);
  sizes[1].ClearSizes();
  sizes[1].AddRepeatedSizes(medSize,1,1);

  InputNode *Win = new InputNode("W input",  sizes, "W", 4);
  

  Sizes ones[2];

  for (Dim dim = 0; dim < 2; ++dim)
    ones[dim].AddRepeatedSizes(one, 1, 1);

  DistType epDist;
  epDist.SetToScalarNoRep();

  InputNode *epIn = new InputNode("ep input",  ones, epDist, "epsilon", 0);
  //InputNode *epIn = new InputNode("ep input",  ones, "epsilon", 0);

  InputNode *tempIn = new InputNode("Temp input",  sizes, "Accum", 4);

  Contraction *cont1 = new Contraction(DMLAYER,COEFONE,COEFZERO,REAL,"abcd","cdij","abij",(string)"cd");
  cont1->AddInputs(6,
		  Uin,0,
		  T1in,0,
		  tempIn,0);

  Poss *poss1 = new Poss(cont1);
  RealPSet *set1 = new RealPSet(poss1);


  Contraction *cont2 = new Contraction(DMLAYER,COEFONE,COEFONE,REAL,"acik","bcjk","abij",(string)"ck");
  cont2->AddInputs(6,
		   Vin,0,
		   T2in,0,
		   set1->OutTun(0),0);

  Poss *poss2 = new Poss(cont2);
  RealPSet *set2 = new RealPSet(poss2);


  Contraction *cont3 = new Contraction(DMLAYER,COEFONE,COEFONE,REAL,"ijkl","abkl","abij",(string)"kl");
  cont3->AddInputs(6,
		   Win,0,
		   T3in,0,
		   set2->OutTun(0),0);

  Poss *poss3 = new Poss(cont3);
  RealPSet *set3 = new RealPSet(poss3);

  Contraction *cont4 = new Contraction(DMLAYER,COEFONE,COEFZERO,REAL,"abij","abij","", (string)"abij");
  cont4->AddInputs(6,
		   T4in,0,
		   set3->OutTun(0),0,
		   epIn,0);

  Poss *poss4 = new Poss(cont4);
  RealPSet *set4 = new RealPSet(poss4);

  OutputNode *out = new OutputNode("output");
  out->AddInput(set4->OutTun(0),0);

  Poss *outerPoss = new Poss(out,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;

}

RealPSet* MartinsExample2()
{
  Sizes sizes[4];

  //a-d = medium
  //i-l = big

  for (Dim dim = 0; dim < 4; ++dim)
    sizes[dim].AddRepeatedSizes(medSize, 1, 1);

  InputNode *Uin = new InputNode("U input",  sizes, "U", 4);

  sizes[2].ClearSizes();
  sizes[2].AddRepeatedSizes(bigSize,1,1);
  sizes[3].ClearSizes();
  sizes[3].AddRepeatedSizes(bigSize,1,1);
  
  InputNode *Vin = new InputNode("V input",  sizes, "V", 4);
  InputNode *Tin = new InputNode("T input",  sizes, "T", 4);


  sizes[0].ClearSizes();
  sizes[0].AddRepeatedSizes(medSize,1,1);
  sizes[1].ClearSizes();
  sizes[1].AddRepeatedSizes(medSize,1,1);

  InputNode *Win = new InputNode("W input",  sizes, "W", 4);
  

  Sizes ones[2];

  for (Dim dim = 0; dim < 2; ++dim)
    ones[dim].AddRepeatedSizes(one, 1, 1);

  DistType epDist;
  epDist.SetToScalarNoRep();

  InputNode *epIn = new InputNode("ep input",  ones, epDist, "epsilon", 0);
  //InputNode *epIn = new InputNode("ep input",  ones, "epsilon", 0);

  InputNode *tempIn = new InputNode("Temp input",  sizes, "Accum", 4);

  Contraction *cont1 = new Contraction(DMLAYER,COEFONE,COEFZERO,REAL,"abcd","cdij","abij",(string)"cd");
  cont1->AddInputs(6,
		  Uin,0,
		  Tin,0,
		  tempIn,0);


  Poss *poss1 = new Poss(cont1);
  RealPSet *set1 = new RealPSet(poss1);

  Contraction *cont2 = new Contraction(DMLAYER,COEFONE,COEFONE,REAL,"acik","bcjk","abij",(string)"ck");
  cont2->AddInputs(6,
		   Vin,0,
		   Tin,0,
		   set1->OutTun(0),0);

  Poss *poss2 = new Poss(cont2);
  RealPSet *set2 = new RealPSet(poss2);


  Contraction *cont3 = new Contraction(DMLAYER,COEFONE,COEFONE,REAL,"ijkl","abkl","abij",(string)"kl");
  cont3->AddInputs(6,
		   Win,0,
		   Tin,0,
		   set2->OutTun(0),0);

  Poss *poss3 = new Poss(cont3);
  RealPSet *set3 = new RealPSet(poss3);

  Contraction *cont4 = new Contraction(DMLAYER,COEFONE,COEFZERO,REAL,"abij","abij","", (string)"abij");
  cont4->AddInputs(6,
		   Tin,0,
		   set3->OutTun(0),0,
		   epIn,0);

  Poss *innerPoss = new Poss(cont4);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *out = new OutputNode("output");
  out->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(out,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;

}

#endif //DOTENSORS
