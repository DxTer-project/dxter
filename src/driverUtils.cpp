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

#include "driverUtils.h"

string DateAndTimeString() {
  time_t t = time(0);
  struct tm* now = localtime(&t);
  return std::to_string((long long int) (now->tm_year + 1900)) + "-"
    + std::to_string((long long int) now->tm_mon + 1) + "-"
    + std::to_string((long long int) now->tm_mday) + "-"
    + std::to_string((long long int) now->tm_hour) + "-"
    + std::to_string((long long int) now->tm_min) + "-"
    + std::to_string((long long int) now->tm_sec);
}

Trans CharToTrans(char c) 
{
  switch (c) {
  case('N'):
    return NORMAL;
  case('T'):
    return TRANS;
  case ('C'):
    return CONJTRANS;
  default:
    throw;
  }
}

Tri CharToTri(char c)
{
  switch (c) {
  case('L'):
    return LOWER;
  case('U'):
    return UPPER;
  default:
    throw;
  }
}

Side CharToSide(char c)
{
  switch (c) {
  case('L'):
    return LEFT;
  case('R'):
    return RIGHT;
  default:
    throw;
  }
}

#if DOLLDLA

VecType CharToVecType(char c) {
  switch(c) {
  case('C'):
    return COLVECTOR;
  case('R'):
    return ROWVECTOR;
  default:
    throw;
  }
}

double BestFlopsPerCycle(Type type, ImplementationRuntimeMap &impTimes, double flopCost) {
  double peakFlopsPerCycle = arch->FlopsPerCycle(type);
  double bestFlopsPerCycle = 0;
  ImplementationRuntimeMapIter mit;
  for (mit = impTimes.begin(); mit != impTimes.end(); ++mit) {
    TimeVecIter vit;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      double totalTimeInCycles = *vit;
      double actualFlopsPerCycle = flopCost / totalTimeInCycles;
      double pctPeak = (actualFlopsPerCycle / peakFlopsPerCycle) * 100;
      if (actualFlopsPerCycle > bestFlopsPerCycle) {
	bestFlopsPerCycle = actualFlopsPerCycle;
      }
      if (pctPeak > 100) {
	cout << "pctPeak > 100\n";
	throw;
      }
    }
  }
  return bestFlopsPerCycle;
}

GraphNum PrintImpMapStats(Type type, ImplementationRuntimeMap &impTimes, double flopCost) {

  double peakFlopsPerCycle = arch->FlopsPerCycle(type);
  GraphNum bestImpNum = 0;
  double overallBestAvgFlopsPerCycle = 0;

  for (auto mit : impTimes) {
    double avgFlopsPerCycle = 0.0;
    double numRuns = 0.0;

    cout << "IMPLEMENTATION # " << std::to_string((long long int) mit.first) << endl;
    for (auto vit : mit.second) {
      double timeInCycles = vit;
      double actualFlopsPerCycle = flopCost / timeInCycles;
      avgFlopsPerCycle += actualFlopsPerCycle;
      numRuns += 1.0;
    }

    avgFlopsPerCycle = avgFlopsPerCycle / numRuns;

    if (avgFlopsPerCycle > overallBestAvgFlopsPerCycle) {
      overallBestAvgFlopsPerCycle = avgFlopsPerCycle;
      bestImpNum = mit.first;
    }

    cout << "Avg. Flops per cycle = " << std::to_string((long double) avgFlopsPerCycle);
    double avgPctPeak = (avgFlopsPerCycle / peakFlopsPerCycle) * 100;
    cout << "\t%Peak = " << std::to_string((long double) avgPctPeak) << endl;
  }
  cout << "Best avg. flops/cycle achieved: " << std::to_string((long double) overallBestAvgFlopsPerCycle) << endl;
  cout << "Best avg. percent of peak: " << std::to_string((long double) (overallBestAvgFlopsPerCycle / peakFlopsPerCycle) * 100) << endl;
  return bestImpNum;
}

#endif // DOLLDLA

Type CharToType(char c)
{
  switch(c) {
#if DOLLDLA
 case('F'):
    return REAL_SINGLE;
 case('D'):
    return REAL_DOUBLE;
#else
  case('R'):
    return REAL;
#endif // DOLLDLA
  case('C'):
    return COMPLEX;
  default:
    throw;
  }
}
