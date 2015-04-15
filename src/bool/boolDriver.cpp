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



#include "base.h"
#include "costs.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "transform.h"
#include "logging.h"
#include <time.h>
#include <iomanip>
#include <chrono>


#if DOBOOL

#include "bool.h"

RealPSet* Example1();

typedef std::chrono::time_point<std::chrono::system_clock> AccurateTime;

double difftime(AccurateTime &end, AccurateTime &start)
{
  ///  std::chrono::duration<double> elapsed_seconds = end-start;
  return (std::chrono::duration_cast<std::chrono::milliseconds>(end-start)).count()/1000.0;
}


void AddTrans()
{

  //  Universe::AddTrans(RedistNode::GetClass(), new CombineAllGathers, BOOLPHASE);
}

void Usage()
{
  cout << "./driver \n";
  cout <<" arg1 == 0  -> Load from file arg1\n";
  cout <<"         1  -> example1\n";
}

int main(int argc, const char* argv[])
{
#ifdef _OPENMP
  omp_set_nested(true);
  omp_init_lock(&RealPSet::m_lock);
#endif
  LOG_START("tensors");
  //  PrintType printType = CODE;
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
      algFunc = Example1;
      break;
    default:
      Usage();
      return 0;
    }
  }

  AddTrans();

  Universe uni;
  AccurateTime start, start2, end;
  uni.PrintStats();

  if (algNum==0) {
    throw;
  }
  else {
    RealPSet *startSet = algFunc();
    uni.Init(startSet);
    uni.Prop();
    GraphIter graphIter(startSet->m_posses.begin()->second);
    Cost flopCost = graphIter.EvalAndSetBest();
    cout << "*****FLOPS = " << setprecision(15) << flopCost << endl;
    start = std::chrono::system_clock::now();
  }


#if BOOlPHASE
    start2 = std::chrono::system_clock::now();
    cout << "Expanding phase 1\n";
    uni.Expand(-1, BOOLPHASE, NULL);
    end = std::chrono::system_clock::now();
    cout << "First phase took " << difftime(end,start2) << " seconds\n";

    cout << "Propagating\n";
    cout.flush();
    start2 = std::chrono::system_clock::now();
    uni.Prop();
    end = std::chrono::system_clock::now();
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif


  end = std::chrono::system_clock::now();
  cout << "Full expansion took " << difftime(end,start) << " seconds\n";
  cout.flush();
  
#if 0
  uni.PrintAll(algNum);
#else
  uni.PrintBest();
#endif


  LOG_END();
  return 0;
}

RealPSet* Example1()
{
  True *tVal = new True;
  False *fVal = new False;
  Or *orNode = new Or;
  orNode->AddInputs0(2, tVal, fVal);
  
  Poss *poss = new Poss(1, orNode);
  RealPSet *pset = new RealPSet(poss);
  return pset;
}


#endif


