#include "row_stride_lldla_primitives.h"
#include "utils.h"

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

  int m = 4;
  int n = m;
  int k = m;
  
  int a_rs = 6;
  int b_rs = 8;
  int c_rs = 4;

  // Do operations and compare results
  dxt_gemm(m, n, k, a_buf, a_rs, b_buf, b_rs, c_buf, c_rs);
  simple_mmul(m, n, k, a_buf, a_rs, 1, b_buf, b_rs, 1, c_buf_copy, c_rs, 1);
  test_buffer_diff(size, c_buf, c_buf_copy, "DXT gemm test");

  return 0;
}
