#include <immintrin.h>
#include "utils.h"
#include <string.h>
#include <unistd.h>

#include "dcaxpy.h"

#define alphaNumRows 1
#define alphaNumCols 1
#define alphaRowStride 1
#define alphaColStride 512
#define xNumRows 512
#define xNumCols 1
#define xRowStride 1
#define xColStride 512
#define yNumRows 512
#define yNumCols 1
#define yRowStride 1
#define yColStride 512
#define MIN_CYCLES 100000000
#define min(a,b) ((a) < (b) ? (a) : (b))

#define MUVALUE 4

typedef union {
  __m256d v;
  double d[4];
} dvec_reg;

void dxt_saxpy_2(double *alpha, double *x, double *y) {
  /*** Algorithm 2 ***
       Unique Num: 2
       Child of: 0
       Result of transformations:
       SVMul register arith - Col vector 0 to 0
       VAddToRegArith 0 to 0
       Cost = 60
  *****************************************/
  dvec_reg alpha_regDup;
  dvec_reg x1_regs;

  dvec_reg y1_regs;

  int lcv0;

  double *x1;
  double *y1;
  //**** (out of 3)
  //------------------------------------//

  alpha_regDup.v = _mm256_broadcast_sd( alpha );
  //**** (out of 1)
  //**** Is real	0 shadows
  x1 = x;
  y1 = y;
  for( lcv0 = xNumRows; lcv0 > 0; lcv0 -= (1*MUVALUE) ) {
    const unsigned int numRowsx = ( (1*MUVALUE) );
    const unsigned int numRowsy = ( (1*MUVALUE) );
    //------------------------------------//

    x1_regs.v = _mm256_load_pd( x1 );
    y1_regs.v = _mm256_load_pd( y1 );
    x1_regs.v = _mm256_mul_pd( alpha_regDup.v , x1_regs.v );
    _mm256_store_pd( x1, x1_regs.v );
    x1_regs.v = _mm256_load_pd( x1 );
    y1_regs.v = _mm256_add_pd( x1_regs.v , y1_regs.v );
    _mm256_store_pd( y1, y1_regs.v );

    //------------------------------------//

    x1 += (1*MUVALUE);
    y1 += (1*MUVALUE);
  }
  //****

  //------------------------------------//

  //****

  /*****************************************/
  // numAlgs = 3
  return;
}
