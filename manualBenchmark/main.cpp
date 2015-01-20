#include <assert.h>
#include <iostream>

#include "dcaxpy.h"
#include "utils.h"

#define NUM_SIZE sizeof(double)
#define ALLOC_BUFFER(size) alloc_aligned_16((size))
#define FILL_WITH_RAND_VALUES(size, buf) rand_doubles((size), (buf))
#define BUF_SIZE 1000000
#define DFLOPS_PER_CYCLE 16
#define SFLOPS_PER_CYCLE 32
#define FLOP_COST 1024
#define MIN_CYCLES 100000000

int main() {
  std::cout << "-------------------- Starting Manual Benchmark --------------------" << std::endl;

  double *alpha = (double*)ALLOC_BUFFER(BUF_SIZE * sizeof(NUM_SIZE));
  double *x = (double*) ALLOC_BUFFER(BUF_SIZE * sizeof(NUM_SIZE));
  double *y = (double*) ALLOC_BUFFER(BUF_SIZE * sizeof(NUM_SIZE));

  FILL_WITH_RAND_VALUES(BUF_SIZE, alpha);
  FILL_WITH_RAND_VALUES(BUF_SIZE, x);
  FILL_WITH_RAND_VALUES(BUF_SIZE, y);


  long long startTime, endTime, totalCycles = 0;
  unsigned int numRuns = 0;
  while (totalCycles < MIN_CYCLES) {
    startTime = rdtsc();
    dxt_saxpy_2(alpha, x, y);
    endTime = rdtsc();
    totalCycles += endTime - startTime;
    numRuns++;
  }

  assert(totalCycles >= MIN_CYCLES);
  
  double avgFlopsPerCycle = (FLOP_COST * numRuns) / ((double) totalCycles);
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
