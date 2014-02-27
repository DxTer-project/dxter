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

#if DOTENSORS

#include "debug.h"


//These control which transformations are included
#define USEVR 1
#define USEMCMR 1
#define USESTAR 1
#define USECONTRIB 1
#define USELOCALCOMP 1
#define EXPLOREREDISTS 1
#define EXPLORETRANS 1
#define USESPECIALTRSM 0
#define USELOWERING 1

#define REMOVESCALEBYONE 1 

//good
#define Hetrmm1 1
//bad
#define Hetrmm2 0
#define Hetrmm3 0
//bad
#define TriInv1 0
#define TriInv2 0
//good
#define TriInv3 1
//scalapack
#define TriInv8 0

#if DOELEM
Size smallSize = 500;
Size medSize = 20000;
Size bigSize = 80000;
Size bs = ELEM_BS;
#elif DOSQM || DOSM
Size smallSize = 500;
Size medSize = 8000;
Size bigSize = 10000;
//Size bs = ELEM_BS;
#endif

void AddTrans()
{

}

void AddSimplifiers()
{ 

}

void Usage()
{
}

int main(int argc, const char* argv[])
{
  //    omp_set_num_threads(1);
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


#if DODPPHASE
  if (CurrPhase == DPPHASE) {
    cout << "Expanding DP phase\n";
    uni.Expand(-1, DPPHASE, DLACullDP);
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

#if DOROPHASE
  if (CurrPhase == ROPHASE) {
    cout << "Expanding RO phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, ROPHASE, DLACullRO);
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

#if DOSR1PHASE
  if (CurrPhase == SR1PHASE) {
    cout << "Expanding SR1 phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, SR1PHASE, DLACullSR);
    time(&end);
    cout << "SR1 phase took " << difftime(end,start2) << " seconds\n";
    
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOSR2PHASE
  if (CurrPhase == SR2PHASE) {
    cout << "Expanding SR2 phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, SR2PHASE, DLACullSR);
    time(&end);
    cout << "SR2 phase took " << difftime(end,start2) << " seconds\n";
    
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif


#if DOSR3PHASE
  if (CurrPhase == SR3PHASE) {
    cout << "Expanding SR3 phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, SR3PHASE, DLACullSR);
    time(&end);
    cout << "SR3 phase took " << difftime(end,start2) << " seconds\n";
    
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif


#if DOSOPHASE
  if (CurrPhase == SOPHASE) {
    cout << "Expanding SO phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, SOPHASE, DLACullSR);
    time(&end);
    cout << "SO phase took " << difftime(end,start2) << " seconds\n";
    
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOSMPPHASE
  if (CurrPhase == SMPPHASE) {
    cout << "Shared-memory parallelization phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, SMPPHASE, DLACullSR);
    time(&end);
    cout << "SMP phase took " << difftime(end,start2) << " seconds\n";
    
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

  uni.PrintAll(algNum);

  /*  if (whichGraph <= 0)
      uni.PrintAll();
      else
      uni.Print(cout, CODE, whichGraph); */

  return 0;
}

#endif //DOTENSORS
