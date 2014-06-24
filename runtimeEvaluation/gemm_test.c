#include "row_stride_lldla_primitives.h"
#include "utils.h"
#include <time.h>
#include <string.h>

int main() {
  // Allocate 16 aligned storage for test buffers
  int dim = 50;
  int size = dim*dim;
  double *a_buf = alloc_aligned_16(size*sizeof(double));
  double *b_buf = alloc_aligned_16(size*sizeof(double));
  double *c_buf = alloc_aligned_16(size*sizeof(double));
  double *c_buf_copy = alloc_aligned_16(size*sizeof(double));

  // Create random test data
  rand_doubles(size, a_buf);
  rand_doubles(size, b_buf);
  rand_doubles(size, c_buf);
  copy_buffer(size, c_buf, c_buf_copy);

  int m = 1200;
  int n = m;
  int k = 4;
  
  int a_rs = 6;
  int b_rs = 8;
  int c_rs = 4;

  // Do operations and compare results
  int i;
  clock_t begin, end;
  double exec_time;
  for (i = 0; i < NUM_ITERATIONS; i++) {
    begin = clock();
    dxt_gemm(m, n, k, a_buf, a_rs, b_buf, b_rs, c_buf, c_rs);
    end = clock();
    exec_time = (double) (end - begin);
    char exec_time_str[100];
    sprintf(exec_time_str, "%f\n", exec_time);
    write(1, exec_time_str, strlen(exec_time_str));
  }
  return 0;
}

