#include <assert.h>
#include <iostream>

#include "rdtsc.h"

#define DFLOPS_PER_CYCLE 16
#define SFLOPS_PER_CYCLE 32
#define FLOP_COST 100000000
#define MIN_CYCLES 100000000

int main() {
  std::cout << "-------------------- Starting Manual Benchmark --------------------" << std::endl;

  long long startTime, endTime, totalCycles = 0;
  unsigned int numRuns = 0;
  while (totalCycles < MIN_CYCLES) {
    startTime = rdtsc();
    endTime = rdtsc();
    totalCycles += endTime - startTime;
    numRuns++;
  }

  assert(totalCycles >= MIN_CYCLES);
  
  double avgFlopsPerCycle = (FLOP_COST * numRuns) / totalCycles;
  double avgPercentOfPeak = (avgFlopsPerCycle / DFLOPS_PER_CYCLE) * 100;

  std::cout << "&&&&&&&&&&&&&&&&&&&& Benchmark Results &&&&&&&&&&&&&&&&&&&&&" << std::endl;
  std::cout << "Number of Runs:       " << numRuns << std::endl;
  std::cout << "Max flops / cycle     " << DFLOPS_PER_CYCLE << std::endl;
  std::cout << "Operation flop count: " << FLOP_COST << std::endl;
  std::cout << "Avg. flops / cycle:   " << avgFlopsPerCycle << std::endl;
  std::cout << "Avg. % of peak:       " << avgPercentOfPeak << std::endl;
  std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;


  std::cout << "-------------------- Done With Manual Benchmark --------------------" << std::endl;
  return 0;
}
