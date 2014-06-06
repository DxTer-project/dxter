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

#if DOLLDLA

#include "debug.h"

Size one = 1;
Size smallSize = 10;
Size medSize = 100;
Size bigSize = 1000;
//Size bs = ELEM_BS;

PSet* Example1();

void AddTrans()
{

}

void AddSimplifiers()
{ 

}

void Usage()
{
  cout << "./driver arg1 arg2 arg3 arg4\n";
  cout <<" arg1 == 0  -> Load from file arg1\n";
  cout <<"         1  -> Basic Example\n";
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
      algFunc = Example1;
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


#if DOLLDLAPHASE
  if (CurrPhase == LLDLAPHASE) {
    cout << "Expanding LL DLA phase\n";
    uni.Expand(-1, LLDLAPHASE, LLDLACull);
    time(&end);
    cout << "LLDLA phase took " << difftime(end,start) << " seconds\n";

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

PSet* Example1()
{
  InputNode *Ain = new InputNode("A input",  smallSize, smallSize, "A");
  InputNode *Bin = new InputNode("B input",  smallSize, smallSize, "B");
  InputNode *Cin = new InputNode("C input",  smallSize, smallSize, "C");

  PossTunnel *tunA = new PossTunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  PossTunnel *tunB = new PossTunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  PossTunnel *tunC = new PossTunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  Contraction *cont = new Contraction(DMLAYER,COEFONE,COEFONE,REAL,"acd", "cefd", "aef", (string)"cd");
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



#endif //DOLLDLA
