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
#include "helperNodes.h"
#include <time.h>
#include <iomanip>
#include "DLAReg.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "contraction.h"
#include "tensorRedist.h"
#include "tensorPermute.h"


#if DOTENSORS

#include "debug.h"
#include "yaxppx.h"
#include "zaxpby.h"
#include "ccsd.h"


Size one = 1;
//Size bs = ELEM_BS;

RealPSet* RedistExample();
RealPSet* RedistExample2();
RealPSet* MartinsExample();
RealPSet* MP2();
RealPSet* MP3();
RealPSet* W();
RealPSet* X();
RealPSet* U();
RealPSet* Q();
RealPSet* P();
RealPSet* H();
RealPSet* F();
RealPSet* G();
RealPSet* z();

void AddTrans()
{
#if 1
  MultiTrans *trans = new MultiTrans;
  trans->AddTrans(new DistContToLocalContStatC(DM2LAYER, SMLAYER));
  trans->AddTrans(new DistContToLocalContStatASumScatter(DM2LAYER, SMLAYER));
  trans->AddTrans(new DistContToLocalContStatBSumScatter(DM2LAYER, SMLAYER));
  Universe::AddTrans(Contraction::GetClass(), trans, DPTENSORPHASE);
#else
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatC(DM2LAYER, SMLAYER), DPTENSORPHASE);
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatASumScatter(DM2LAYER, SMLAYER), DPTENSORPHASE);
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatBSumScatter(DM2LAYER, SMLAYER), DPTENSORPHASE);  
#endif

  for (Dim dim = 0; dim < 10; ++dim) {
    Universe::AddTrans(Contraction::GetClass(), new ContractionLoopExp(ABSLAYER, DM1LAYER, dim), DPTENSORPHASE);
    Universe::AddTrans(Contraction::GetClass(), new ContractionLoopExp(DM1LAYER, DM2LAYER, dim), DPTENSORPHASE);
  }

  Universe::AddTrans(Contraction::GetClass(), new ContractionLowerLayer(ABSLAYER, DM2LAYER, TensorBS.GetSize()), DPTENSORPHASE);
  Universe::AddTrans(Contraction::GetClass(), new ContractionLowerLayer(DM1LAYER, DM2LAYER, TensorBS.GetSize()), DPTENSORPHASE);

#if 0
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatAAllReduce(DMLAYER, SMLAYER), DPTENSORPHASE);
  Universe::AddTrans(Contraction::GetClass(), new DistContToLocalContStatBAllReduce(DMLAYER, SMLAYER), DPTENSORPHASE);
#endif

  Universe::AddTrans(SumScatterUpdateNode::GetClass(), new SeparateRedistFromSumScatter, SUMSCATTERTENSORPHASE);
  Universe::AddTrans(SumScatterUpdateNode::GetClass(), new MoveSumScatterRedistAfter, SUMSCATTERTENSORPHASE);

  Universe::AddTrans(YAxpPx::GetClass(), new DistYAxpPxToDefaultLocalYAxpPx, DPTENSORPHASE);

  Universe::AddTrans(YAxpPx::GetClass(), new YAxpPxLowerLayer(ABSLAYER,DM1LAYER,TensorBS.GetSize()), DPTENSORPHASE);
  Universe::AddTrans(YAxpPx::GetClass(), new YAxpPxLowerLayer(DM1LAYER,DM2LAYER,TensorBS.GetSize()), DPTENSORPHASE);
  
  Universe::AddTrans(ZAxpBy::GetClass(), new ZAxpByLowerLayer(ABSLAYER,SMLAYER), DPTENSORPHASE);

#if ALLMULTIMODEALLGATHER
    Universe::AddTrans(RedistNode::GetClass(), new SplitAllAllGathers, ROTENSORPHASE);
#endif

#if DOPACKOPTPHASE
    Universe::AddTrans(Contraction::GetClass(), new PermuteWhileUnpacking(0), PACKOPTPHASE);
    Universe::AddTrans(Contraction::GetClass(), new PermuteWhileUnpacking(1), PACKOPTPHASE);
    Universe::AddTrans(Contraction::GetClass(), new PermuteWhileUnpacking(2), PACKOPTPHASE);
#endif

#if DOFINALOPTPHASE
    Universe::AddTrans(LoopTunnel::GetClass(), new PermuteLoopHoist, SIMP);
#endif
  
#if 1
  for(Dim dim = 0; dim < NUM_GRID_DIMS; ++dim) {
    Universe::AddTrans(YAxpPx::GetClass(), new YAxpPxLoopExp(ABSLAYER,DM1LAYER,dim), DPTENSORPHASE);
    Universe::AddTrans(YAxpPx::GetClass(), new YAxpPxLoopExp(DM1LAYER,DM2LAYER,dim), DPTENSORPHASE);
    Universe::AddTrans(RedistNode::GetClass(), new SplitRedistribs(dim), ROTENSORPHASE);
    Universe::AddTrans(RedistNode::GetClass(), new SingleIndexAllToAll(dim), ROTENSORPHASE);
    Universe::AddTrans(RedistNode::GetClass(), new DoubleIndexAllToAll(dim), ROTENSORPHASE);
    Universe::AddTrans(RedistNode::GetClass(), new DoubleIndexAllToAllPrefix(dim), ROTENSORPHASE);
    Universe::AddTrans(RedistNode::GetClass(), new SplitAllGathers(dim), ROTENSORPHASE);

#if SPLITMULTIMODESCATTER
    Universe::AddTrans(SumScatterUpdateNode::GetClass(), new SplitSumScatter(dim), SUMSCATTERTENSORPHASE);
#endif

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
   Universe::AddTrans(Permute::GetClass(), new CombinePermutations, SIMP);
   Universe::AddTrans(Permute::GetClass(), new MovePermuteIntoTempVarNode, SIMP);
   Universe::AddTrans(ScaleNode::GetClass(), new RemoveScaleByOne, SIMP);
   Universe::AddTrans(TempVarNode::GetClass(), new MoveTempVarNodeIntoLoop, SIMP);
}

void Usage()
{
  cout << "./driver arg1 arg2 arg3 arg4\n";
  cout <<" arg1 == 0  -> Load from file arg1\n";
  cout <<"         1  -> Redist example\n";
  cout <<"         2  -> Redist Example2\n";
  cout <<"         3  -> Martin's Example\n";
  cout <<"         4  -> MP2\n";
  cout <<"         5  -> MP3\n";
  cout <<"         6  -> W_bmje\n";
  cout <<"         7  -> X_bmej\n";
  cout <<"         8  -> U_mnie\n";
  cout <<"         9  -> Q_mnij\n";
  cout <<"        10  -> P_jimb\n";
  cout <<"        11  -> H_me\n";
  cout <<"        12  -> F_ae\n";
  cout <<"        13  -> G_mi\n";
  cout <<"        14  -> z_ai\n";
}

int main(int argc, const char* argv[])
{
#ifdef _OPENMP
  omp_set_num_threads(1);
  omp_set_nested(true);
  omp_init_lock(&RealPSet::m_lock);
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
      algFunc = RedistExample;
      break;
    case(2):
      algFunc = RedistExample2;
      break;
    case(3):
      algFunc = MartinsExample;
      break;
    case(4):
      algFunc = MP2;
      break;
    case(5):
      algFunc = MP3;
      break;
    case(6):
      algFunc = W;
      break;
    case(7):
      algFunc = X;
      break;
    case(8):
      algFunc = U;
      break;
    case(9):
      algFunc = Q;
      break;
    case(10):
      algFunc = P;
      break;
    case(11):
      algFunc = H;
      break;
    case(12):
      algFunc = F;
      break;
    case(13):
      algFunc = G;
      break;
    case(14):
      algFunc = z;
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
    RealPSet *startSet = algFunc();
    uni.Init(startSet);
    uni.Prop();
    GraphIter graphIter(startSet->m_posses.begin()->second);
    //    cout << "Printing evaluation code\n";
    Cost flopCost = graphIter.EvalAndSetBest();
    cout << "*****FLOPS = " << setprecision(15) << flopCost << endl;
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
    uni.CullWorstPerformers(.99, 3);
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOFUSEANDOPTPHASE
  if (CurrPhase == FUSEANDOPTTENSORPHASE) {
    cout << "Fusing and opt phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, FUSEANDOPTTENSORPHASE, TenCullRO);
    time(&end);
    cout << "Fusing and opt phase took " << difftime(end,start2) << " seconds\n";
    
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif


#if DOPACKOPTPHASE
  if (CurrPhase == PACKOPTPHASE) {
    cout << "Pack optimization phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Prop();
    uni.CullWorstPerformers(.99, 3);
    time(&end);
    cout << "After culling worst (" << difftime(end,start2) << " secs), left with " << uni.TotalCount() << endl;
    time(&start2);
    uni.InlineAllSets();
    time(&end);
    cout << "Inlining took " << difftime(end,start2) << " seconds\n";
    time(&start2);
    uni.Expand(numIters, PACKOPTPHASE, TenCullRO);
    time(&end);
    cout << "Pack optimization phase took " << difftime(end,start2) << " seconds\n";

    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOFINALOPTPHASE
  if (CurrPhase == FINALOPTPHASE) {
    cout << "Final optimization phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Prop();
    uni.CullAllBut(1);
    time(&end);
    cout << "After culling worst (" << difftime(end,start2) << " secs), left with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, FINALOPTPHASE, TenCullRO);
    time(&end);
    cout << "Pack optimization phase took " << difftime(end,start2) << " seconds\n";

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

RealPSet* RedistExample()
{
  Size bigSize = 1000;

  Sizes sizes[4];

  for (Dim dim = 0; dim < 4; ++dim)
    sizes[dim].AddRepeatedSizes(bigSize, 1, 1);

  InputNode *Ain = new InputNode("A input",  sizes, "A", 4);

  DistType type1;
  type1.SetToDefault(4);
  type1.m_dists[0].m_val = 3;
  type1.m_dists[1].m_val = 4;
  type1.m_dists[2].m_val = 0;
  type1.m_dists[3].m_val = 0;

  DimVec ident;
  IdentDimVec(4, ident);
  
  RedistNode *redist1 = new RedistNode(type1, Ain->GetNameStr(0), ident, ident);
  redist1->AddInput(Ain, 0);

  DistType type2;
  type2.SetToDefault(4);
  type2.m_dists[0].m_val = 0;
  type2.m_dists[1].m_val = 2;
  type2.m_dists[2].m_val = 0;
  type2.m_dists[3].m_val = 4;

  RedistNode *redist2 = new RedistNode(type2, Ain->GetNameStr(0), ident, ident);
  redist2->AddInput(Ain, 0);

  OutputNode *Cout1 = new OutputNode("C output");
  Cout1->AddInput(redist1, 0);


  OutputNode *Cout2 = new OutputNode("C1 output");
  Cout2->AddInput(redist2, 0);

  Poss *outerPoss = new Poss(2, Cout1, Cout2);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;

}

RealPSet* RedistExample2()
{
  Size bigSize = 1000;
  
  Sizes sizes[4];

  for (Dim dim = 0; dim < 4; ++dim)
    sizes[dim].AddRepeatedSizes(bigSize, 1, 1);

  InputNode *Ain = new InputNode("A input",  sizes, "A", 4);

  DistType type1;
  type1.SetToDefault(4);
  type1.m_dists[0].m_val = 0;
  type1.m_dists[1].m_val = 2;
  type1.m_dists[2].m_val = 0;
  type1.m_dists[3].m_val = 4;

  DimVec ident;
  IdentDimVec(4, ident);

  RedistNode *redist1 = new RedistNode(type1, Ain->GetNameStr(0), ident, ident);
  redist1->AddInput(Ain, 0);

  OutputNode *Cout1 = new OutputNode("C output");
  Cout1->AddInput(redist1, 0);

  Poss *outerPoss = new Poss(1, Cout1);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}


RealPSet* MartinsExample()
{
  Size medSize = 100;
  Size bigSize = 1000;

  Sizes sizes[4];

  //a-d = medium
  //i-l = big

  for (Dim dim = 0; dim < 4; ++dim)
    sizes[dim].AddRepeatedSizes(bigSize, 1, 1);

  InputNode *Uin = new InputNode("U input",  sizes, "U", 4);
  

  sizes[2].ClearSizes();
  sizes[2].AddRepeatedSizes(medSize,1,1);
  sizes[3].ClearSizes();
  sizes[3].AddRepeatedSizes(medSize,1,1);
  
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

  sizes[0].ClearSizes();
  sizes[0].AddRepeatedSizes(bigSize,1,1);
  sizes[1].ClearSizes();
  sizes[1].AddRepeatedSizes(bigSize,1,1);

  InputNode *tempIn = new InputNode("Temp input",  sizes, "Accum", 4);

  Contraction *cont1 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"abcd","cdij","abij",(string)"cd");
  cont1->AddInputs(6,
		  Uin,0,
		  Tin,0,
		  tempIn,0);


  Poss *poss1 = new Poss(cont1);
  RealPSet *set1 = new RealPSet(poss1);

  Contraction *cont2 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"acik","bcjk","abij",(string)"ck");
  cont2->AddInputs(6,
		   Vin,0,
		   Tin,0,
		   set1->OutTun(0),0);

  Poss *poss2 = new Poss(cont2);
  RealPSet *set2 = new RealPSet(poss2);


  Contraction *cont3 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"ijkl","abkl","abij",(string)"kl");
  cont3->AddInputs(6,
		   Win,0,
		   Tin,0,
		   set2->OutTun(0),0);

  Poss *poss3 = new Poss(cont3);
  RealPSet *set3 = new RealPSet(poss3);

  Contraction *cont4 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"abij","abij","", (string)"abij");
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

RealPSet* MP2()
{
  Size eSize = 53*5;
  Size fSize = 53*5;
  Size mSize = 5*5;
  Size nSize = 5*5;


  InputNode *v_efmn;
  InputNode *t_efmn;


  InputNode *axppx1_temp;

  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(eSize,1,1);
    sizes[1].AddRepeatedSizes(fSize,1,1);
    sizes[2].AddRepeatedSizes(mSize,1,1);
    sizes[3].AddRepeatedSizes(nSize,1,1);
    t_efmn = new InputNode("t_efmn", sizes, "t_efmn", 4);
  }

  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(eSize,1,1);
    sizes[1].AddRepeatedSizes(fSize,1,1);
    sizes[2].AddRepeatedSizes(mSize,1,1);
    sizes[3].AddRepeatedSizes(nSize,1,1);
    v_efmn = new InputNode("v_efmn", sizes, "v_efmn", 4);
  }

  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(eSize,1,1);
    sizes[1].AddRepeatedSizes(fSize,1,1);
    sizes[2].AddRepeatedSizes(mSize,1,1);
    sizes[3].AddRepeatedSizes(nSize,1,1);
    axppx1_temp = new InputNode("axppx1_temp", 
				sizes, "axppx1_temp", 4);
  }


  Sizes ones[2];
  for (Dim dim = 0; dim < 2; ++dim)
    ones[dim].AddRepeatedSizes(one, 1, 1);

  DistType scalarDist;
  scalarDist.SetToScalarNoRep();

  InputNode *scalarIn = new InputNode("scalar input",  ones, scalarDist, "E_MP2", 0);


  YAxpPx *axppx1 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONE, "efnm", "efmn");
  axppx1->AddInputs(6,
		   v_efmn, 0,
		   v_efmn, 0,
		   axppx1_temp, 0);
  Poss *axppx1Poss = new Poss(axppx1);
  RealPSet * axppx1Set = new RealPSet(axppx1Poss);


  Contraction *cont1 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"efmn","efmn","",(string)"efmn");
  cont1->AddInputs(6,
		   axppx1Set->OutTun(0),0,
		   t_efmn,0,
		   scalarIn,0);
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  OutputNode *out = new OutputNode("output");
  out->AddInput(cont1Set->OutTun(0),0);

  Poss *outerPoss = new Poss(out, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MP3()
{
  const Size bigMP3Size = 300;
  const Size smallMP3Size = 300;
  Size eSize = bigMP3Size;
  Size fSize = bigMP3Size;
  Size gSize = bigMP3Size;
  Size hSize = bigMP3Size;
  Size mSize = smallMP3Size;
  Size nSize = smallMP3Size;
  Size oSize = smallMP3Size;
  Size pSize = smallMP3Size;


  InputNode *t_efmn;
  InputNode *v_opmn;
  InputNode *v_efgh;
  InputNode *v_oegm;
  InputNode *v2_oegm;
  InputNode *accum_temp;
  InputNode *cont1_temp;
  InputNode *axppx2_temp;
  InputNode *axppx3_temp;

  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(eSize,1,1);
    sizes[1].AddRepeatedSizes(fSize,1,1);
    sizes[2].AddRepeatedSizes(mSize,1,1);
    sizes[3].AddRepeatedSizes(nSize,1,1);
    t_efmn = new InputNode("t_efmn", sizes, "t_efmn", 4);
  }

  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(oSize,1,1);
    sizes[1].AddRepeatedSizes(pSize,1,1);
    sizes[2].AddRepeatedSizes(mSize,1,1);
    sizes[3].AddRepeatedSizes(nSize,1,1);
    v_opmn = new InputNode("v_opmn", sizes, "v_opmn", 4);
  }


  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(eSize,1,1);
    sizes[1].AddRepeatedSizes(fSize,1,1);
    sizes[2].AddRepeatedSizes(gSize,1,1);
    sizes[3].AddRepeatedSizes(hSize,1,1);
    v_efgh = new InputNode("v_efgh", sizes, "v_efgh", 4);
  }



  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(oSize,1,1);
    sizes[1].AddRepeatedSizes(eSize,1,1);
    sizes[2].AddRepeatedSizes(gSize,1,1);
    sizes[3].AddRepeatedSizes(mSize,1,1);
    v_oegm = new InputNode("v_oegm", sizes, "v_oegm", 4);
  }

  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(oSize,1,1);
    sizes[1].AddRepeatedSizes(eSize,1,1);
    sizes[2].AddRepeatedSizes(gSize,1,1);
    sizes[3].AddRepeatedSizes(mSize,1,1);
    v2_oegm = new InputNode("v2_oegm", sizes, "v2_oegm", 4);
  }



  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(eSize,1,1);
    sizes[1].AddRepeatedSizes(fSize,1,1);
    sizes[2].AddRepeatedSizes(mSize,1,1);
    sizes[3].AddRepeatedSizes(nSize,1,1);
    accum_temp = new InputNode("accum_temp", 
			       sizes, "accum_temp", 4);
  }

  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(eSize,1,1);
    sizes[1].AddRepeatedSizes(fSize,1,1);
    sizes[2].AddRepeatedSizes(mSize,1,1);
    sizes[3].AddRepeatedSizes(nSize,1,1);
    cont1_temp = new InputNode("cont1_temp", 
			       sizes, "cont1_temp", 4);
  }

  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(gSize,1,1);
    sizes[1].AddRepeatedSizes(fSize,1,1);
    sizes[2].AddRepeatedSizes(oSize,1,1);
    sizes[3].AddRepeatedSizes(nSize,1,1);
    axppx2_temp = new InputNode("axppx2_temp", 
			       sizes, "axppx2_temp", 4);
  }

  {
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(oSize,1,1);
    sizes[1].AddRepeatedSizes(eSize,1,1);
    sizes[2].AddRepeatedSizes(gSize,1,1);
    sizes[3].AddRepeatedSizes(mSize,1,1);
    axppx3_temp = new InputNode("axppx3_temp", 
			       sizes, "axppx3_temp", 4);
  }






  Sizes ones[2];
  for (Dim dim = 0; dim < 2; ++dim)
    ones[dim].AddRepeatedSizes(one, 1, 1);

  DistType scalarDist;
  scalarDist.SetToScalarNoRep();

  InputNode *scalarIn = new InputNode("scalar input",  ones, scalarDist, "E_MP3", 0);


  Contraction *cont1 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"oegm","gfno","efmn",(string)"go");
  cont1->AddInputs(6,
		  v2_oegm,0,
		  t_efmn,0,
		  cont1_temp,0);
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  YAxpPx *axppx1 = new YAxpPx(ABSLAYER, COEFONEHALF, COEFONE, "efmn", "efnm");
  axppx1->AddInputs(6,
		   cont1Set->OutTun(0), 0,
		   cont1Set->OutTun(0), 0,
		   accum_temp, 0);
  Poss *axppx1Poss = new Poss(axppx1);
  RealPSet * axppx1Set = new RealPSet(axppx1Poss);

  YAxpPx *axppx2 = new YAxpPx(ABSLAYER,COEFTWO, COEFNEGONE, "gfon", "gfno");
  axppx2->AddInputs(6,
		   t_efmn, 0,
		   t_efmn, 0,
		   axppx2_temp, 0);
  Poss *axppx2Poss = new Poss(axppx2);
  RealPSet * axppx2Set = new RealPSet(axppx2Poss);

  ZAxpBy *axppx3 = new ZAxpBy(ABSLAYER,COEFTWO, COEFNEGONE);
  axppx3->AddInputs(6,
		   v_oegm, 0,
		   v2_oegm, 0,
		   axppx3_temp, 0);
  Poss *axppx3Poss = new Poss(axppx3);
  RealPSet * axppx3Set = new RealPSet(axppx3Poss);

  Contraction *cont2 = new Contraction(ABSLAYER,COEFONEHALF,COEFNEGONE,REAL,"oegm","gfon","efmn",(string)"go");
  cont2->AddInputs(6,
		   axppx3Set->OutTun(0),0,
		   axppx2Set->OutTun(0),0,
		   axppx1Set->OutTun(0),0);
  Poss *cont2Poss = new Poss(cont2);
  RealPSet *cont2Set = new RealPSet(cont2Poss);

  Contraction *cont3 = new Contraction(ABSLAYER,COEFONEHALF,COEFONE,REAL,"efgh","ghmn","efmn",(string)"gh");
  cont3->AddInputs(6,
		   v_efgh, 0,
		   t_efmn, 0,
		   cont2Set->OutTun(0),0);
  Poss *cont3Poss = new Poss(cont3);
  RealPSet *cont3Set = new RealPSet(cont3Poss);

  Contraction *cont4 = new Contraction(ABSLAYER,COEFONEHALF,COEFONE,REAL,"opmn","efop","efmn",(string)"op");
  cont4->AddInputs(6,
		   v_opmn, 0,
		   t_efmn, 0,
		   cont3Set->OutTun(0),0);
  Poss *cont4Poss = new Poss(cont4);
  RealPSet *cont4Set = new RealPSet(cont4Poss);


  Contraction *cont5 = new Contraction(ABSLAYER,COEFTWO,COEFZERO,REAL,"efmn","efmn","",(string)"efmn");
  cont5->AddInputs(6,
		   axppx2Set->OutTun(0), 0,
		   cont4Set->OutTun(0), 0,
		   scalarIn,0);
  Poss *cont5Poss = new Poss(cont5);
  RealPSet *cont5Set = new RealPSet(cont5Poss);


  OutputNode *out = new OutputNode("output");
  out->AddInput(cont5Set->OutTun(0),0);

  Poss *outerPoss = new Poss(out, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* W()
{
  //~ 10:1 ratio
  // 53, 5 for H20
  const Size big = 530; //a-h
  const Size small = 50; //i-p

  InputNode *w_bmje = CreateInput4("w_bmje", big, small, small, big);
  InputNode *x_bmej = CreateInput4("x_bmej", big, small, big, small);
  InputNode *r_bmef = CreateInput4("r_bmfe", big, small, big, big);
  InputNode *t_fj = CreateInput2("t_fj", big, small);
  InputNode *u_mnje = CreateInput4("u_mnje", small, small, small, big);
  InputNode *v_femn = CreateInput4("v_femn", big, big, small, small);
  InputNode *T_bfnj = CreateInput4("T_bfnj", big, big, small, small);

  RealPSet *set = W_bmje_calc(w_bmje, x_bmej, r_bmef, t_fj, 
			      u_mnje, v_femn, T_bfnj, 
			      big, small);

  OutputNode *out = new OutputNode("output");
  out->AddInput(set->OutTun(0),0);

  Poss *outerPoss = new Poss(out, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* X()
{
  //~ 10:1 ratio
  // 530, 5 for H20
  const Size big = 530; //a-h
  const Size small = 50; //i-p

  InputNode *x_bmej = CreateInput4("x_bmej", big, small, big, small);
  InputNode *r_bmef = CreateInput4("r_bmfe", big, small, big, big);
  InputNode *t_fj = CreateInput2("t_fj", big, small);
  InputNode *u_mnje = CreateInput4("u_mnje", small, small, small, big);
  InputNode *v_femn = CreateInput4("v_femn", big, big, small, small);
  InputNode *T_bfnj = CreateInput4("T_bfnj", big, big, small, small);

  RealPSet *set = X_bmej_calc(x_bmej, r_bmef, t_fj, 
			      u_mnje, v_femn, T_bfnj, 
			      big, small);

  OutputNode *out = new OutputNode("output");
  out->AddInput(set->OutTun(0),0);

  Poss *outerPoss = new Poss(out, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* U()
{
  //~ 10:1 ratio
  // 530, 5 for H20
  const Size big = 530; //a-h
  const Size small = 50; //i-p

  InputNode *t_fj = CreateInput2("t_fj", big, small);
  InputNode *u_mnie = CreateInput4("u_mnie", small, small, small, big);
  InputNode *v_femn = CreateInput4("v_femn", big, big, small, small);

  RealPSet *set = U_mnie_calc(t_fj, 
			      u_mnie, v_femn,
			      big, small);

  OutputNode *out = new OutputNode("output");
  out->AddInput(set->OutTun(0),0);

  Poss *outerPoss = new Poss(out, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* Q()
{
  //~ 10:1 ratio
  // 530, 5 for H20
  const Size big = 530; //a-h
  const Size small = 50; //i-p

  InputNode *q_mnij = CreateInput4("q_mnij", small, small, small, small);
  InputNode *t_fj = CreateInput2("t_fj", big, small);
  InputNode *u_mnie = CreateInput4("u_mnie", small, small, small, big);
  InputNode *v_femn = CreateInput4("v_femn", big, big, small, small);
  InputNode *T_bfnj = CreateInput4("T_bfnj", big, big, small, small);

  RealPSet *set = Q_mnij_calc(q_mnij, t_fj, 
			      u_mnie, v_femn, T_bfnj, 
			      big, small);

  OutputNode *out = new OutputNode("output");
  out->AddInput(set->OutTun(0),0);

  Poss *outerPoss = new Poss(out, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* P()
{
  //~ 10:1 ratio
  // 530, 5 for H20
  const Size big = 530; //a-h
  const Size small = 50; //i-p

  InputNode *u_jimb = CreateInput4("u_jimb", small, small, small, big);
  InputNode *r_bmef = CreateInput4("r_bmef", big, small, big, big);
  InputNode *t_fj = CreateInput2("t_fj", big, small);
  InputNode *w_bmie = CreateInput4("w_bmie", big, small, small, big);
  InputNode *T_efij = CreateInput4("T_efij", big, big, small, small);
  InputNode *x_bmej = CreateInput4("x_bmej", big, small, big, small);

  RealPSet *set = P_jimb_calc(r_bmef, t_fj, 
			      u_jimb, w_bmie, T_efij, 
			      x_bmej,
			      big, small);

  OutputNode *out = new OutputNode("output");
  out->AddInput(set->OutTun(0),0);

  Poss *outerPoss = new Poss(out, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* H()
{
  //~ 10:1 ratio
  // 530, 5 for H20
  const Size big = 530; //a-h
  const Size small = 50; //i-p

  InputNode *t_fn = CreateInput2("t_fn", big, small);
  InputNode *v_efmn = CreateInput4("v_efmn", big, big, small, small);

  RealPSet *set = H_me_calc(t_fn, v_efmn,
			      big, small);

  OutputNode *out = new OutputNode("output");
  out->AddInput(set->OutTun(0),0);

  Poss *outerPoss = new Poss(out, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}


RealPSet* F()
{
  //~ 10:1 ratio
  // 530, 5 for H20
  const Size big = 530; //a-h
  const Size small = 50; //i-p
  
  InputNode *H_me = CreateInput2("H_me", small, big);
  InputNode *T_afmn = CreateInput4("T_afmn", big, big, small, small);
  InputNode *r_amef = CreateInput4("r_amef", big, small, big, big);
  InputNode *t_am = CreateInput2("t_am", big, small);
  InputNode *v_efmn = CreateInput4("v_efmn", big, big, small, small);
  
  RealPSet *set = F_ae_calc(H_me, r_amef, 
			    t_am, v_efmn, T_afmn,
			    big, small);
  
  OutputNode *out = new OutputNode("output");
  out->AddInput(set->OutTun(0),0);

  Poss *outerPoss = new Poss(out, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}



RealPSet* G()
{
  //~ 10:1 ratio
  // 530, 5 for H20
  const Size big = 530; //a-h
  const Size small = 50; //i-p
  
  InputNode *H_me = CreateInput2("H_me", small, big);
  InputNode *T_efin = CreateInput4("T_efin", big, big, small, small);
  InputNode *u_mnie = CreateInput4("u_mnie", small, small, small, big);
  InputNode *t_ei = CreateInput2("t_ei", big, small);
  InputNode *v_efmn = CreateInput4("v_efmn", big, big, small, small);
  
  RealPSet *set = G_mi_calc(H_me, u_mnie,
			    t_ei, v_efmn, T_efin,
			    big, small);
  
  OutputNode *out = new OutputNode("output");
  out->AddInput(set->OutTun(0),0);

  Poss *outerPoss = new Poss(out, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* z()
{
  //~ 10:1 ratio
  // 530, 5 for H20
  const Size big = 530; //a-h
  const Size small = 50; //i-p
  
  InputNode *f_ae = CreateInput2("f_ae", big, big);
  InputNode *G_mi = CreateInput2("G_mi", small, small);
  InputNode *w_amie = CreateInput4("w_amie", big, small, small, big);
  InputNode *x_amei = CreateInput4("x_amei", big, small, big, small);
  InputNode *T_aeim = CreateInput4("T_aeim", big, big, small, small);
  InputNode *r_amef = CreateInput4("r_amef", big, small, big, big);
  InputNode *t_am = CreateInput2("t_am", big, small);
  DLANode *H_me = CreateInput2("H_me", small, big);
  InputNode *U_mnie = CreateInput4("U_mnie", small, small, small, big);
  

  RealPSet *set = z_ai_calc(f_ae, G_mi, 
			    H_me, U_mnie,
			    w_amie,
			    x_amei, 
			    t_am, r_amef, T_aeim,
			    big, small);
  
  OutputNode *out = new OutputNode("output");
  out->AddInput(set->OutTun(0),0);

  Poss *outerPoss = new Poss(out, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}



#endif //DOTENSORS
