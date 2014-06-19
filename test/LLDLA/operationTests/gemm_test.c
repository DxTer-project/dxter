#include "row_stride_lldla_primitives.h"
#include "utils.h"

#define MUVALUE 2

void dxt_gemm(int CNumCols, int CNumRows, int ANumCols,
	      double *A, int ARowStride,
	      double *B, int BRowStride,
	      double *C, int CRowStride) {
  
  double b = 1.0;
  double *beta = &b;
  double *A1, *B1, *C1, *A11, *B11, *C11;
  //**** (out of 4)
  //------------------------------------//

  B1 = B;
  C1 = C;
  //Dim-n loop
  while ( C1 < (C + CNumCols) )  {
    //**** (out of 2)
    //------------------------------------//

    A1 = A;
    C11 = C1;
    //Dim-m loop
    while ( C11 < (C1 + CNumRows * CRowStride) )  {
      //**** (out of 1)
      //------------------------------------//

      row_stride_smul_2x2( beta, C11, CRowStride);
      A11 = A1;
      B11 = B1;
      //Dim-k loop
      while ( A11 < (A1 + ANumCols) )  {
	//**** (out of 1)
	//------------------------------------//

	row_stride_mmul_2x2_2x2( A11, ARowStride, B11, BRowStride, C11, CRowStride);

	//------------------------------------//

	//****
	A11 += MUVALUE;
	B11 += BRowStride * MUVALUE;
      }

      //------------------------------------//

      //****
      A1 += ARowStride * MUVALUE;
      C11 += CRowStride * MUVALUE;
    }

    //------------------------------------//

    //****
    B1 += MUVALUE;
    C1 += MUVALUE;
  }

  //------------------------------------//

  //****

}

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
