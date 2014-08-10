#ifndef COL_STRIDE_LLDLA_PRIMITIVES_H
#define COL_STRIDE_LLDLA_PRIMITIVES_H

#include "utils.h"

#define G(i, j) a[ (i) + (j)*a_col_stride ]
#define H(i, j) b[ (i) + (j)*b_col_stride ]
#define I(i, j) c[ (i) + (j)*c_col_stride ]

// Addition kernels

void col_stride_add_2x1(const double *a, int a_col_stride,
			double *b, int b_col_stride) {
  v2df_t b_00_10;
  b_00_10.v = VEC_PD_LOAD(&H(0, 0));
  b_00_10.v += VEC_PD_LOAD(&G(0, 0));
  VEC_PD_STORE(&H(0, 0), b_00_10.v);
}

void col_stride_add_1x2(const double *a, int a_col_stride,
			double *b, int b_col_stride) {
  v2df_t a_00_01, b_00_01;
  a_00_01.v = VEC_D_LOAD(&G(0, 0), &G(0, 1));
  b_00_01.v = VEC_D_LOAD(&H(0, 0), &H(0, 1));
  b_00_01.v += a_00_01.v;
  H(0, 0) = b_00_01.d[0];
  H(0, 1) = b_00_01.d[1];
}

void col_stride_add_2x2(const double *a, int a_col_stride,
			double *b, int b_col_stride) {
  col_stride_add_2x1(&G(0, 0), a_col_stride, &H(0, 0), b_col_stride);
  col_stride_add_2x1(&G(0, 1), a_col_stride, &H(0, 1), b_col_stride);
}

// Scalar multiplication kernels

void col_stride_smul_2x1(const double *scalar,
			 double *a, int a_col_stride) {

  v2df_t scal, a_00_10;
  
  scal.v = VEC_DUP_LOAD(scalar);
  a_00_10.v = VEC_PD_LOAD(&G(0, 0));
  a_00_10.v *= scal.v;
  VEC_PD_STORE(&G(0, 0), a_00_10.v);
}

void col_stride_smul_1x2(const double *scalar,
			 double *a, int a_col_stride) {
  v2df_t scal, a_00_01;
  scal.v = VEC_DUP_LOAD(scalar);
  a_00_01.v = VEC_D_LOAD(&G(0, 0), &G(0, 1));
  a_00_01.v = a_00_01.v * scal.v;
  G(0, 0) = a_00_01.d[0];
  G(0, 1) = a_00_01.d[1];
}

void col_stride_smul_2x2(const double *scalar,
			 double *a, int a_col_stride) {
  col_stride_smul_2x1(scalar, &G(0, 0), a_col_stride);
  col_stride_smul_2x1(scalar, &G(0, 1), a_col_stride);
}

// Matrix multiply primitives

void col_stride_mmul_1x2_2x1(const double *a, int a_col_stride,
			     const double *b, int b_col_stride,
			     double *c, int c_col_stride) {
  v2df_t a_00_01, b_00_10;
  
  a_00_01.v = VEC_D_LOAD(&G(0, 0), &G(0, 1));
  b_00_10.v = VEC_PD_LOAD(&H(0, 0));
  a_00_01.v *= b_00_10.v;
  I(0, 0) += a_00_01.d[0];
  I(0, 0) += a_00_01.d[1];
}

void col_stride_mmul_2x1_1x2(const double *a, int a_col_stride,
			     const double *b, int b_col_stride,
			     double *c, int c_col_stride) {

  v2df_t a_00_10, b_00, b_01;

  a_00_10.v = VEC_PD_LOAD(&G(0, 0));
  b_00.v = VEC_DUP_LOAD(&H(0, 0));
  b_01.v = VEC_DUP_LOAD(&H(0, 1));

  b_00.v *= a_00_10.v;
  b_01.v *= a_00_10.v;

  I(0, 0) += b_00.d[0];
  I(1, 0) += b_00.d[1];
  I(0, 1) += b_01.d[0];
  I(1, 1) += b_01.d[1];
}

void col_stride_mmul_1x2_2x2(const double *a, int a_col_stride,
			     const double *b, int b_col_stride,
			     double *c, int c_col_stride) {

  col_stride_mmul_1x2_2x1(&G(0, 0), a_col_stride,
			  &H(0, 0), b_col_stride,
			  &I(0, 0), c_col_stride);

  col_stride_mmul_1x2_2x1(&G(0, 0), a_col_stride,
			  &H(0, 1), b_col_stride,
			  &I(0, 1), c_col_stride);
}

void col_stride_mmul_2x2_2x1(const double *a, int a_col_stride,
			     const double *b, int b_col_stride,
			     double *c, int c_col_stride) {

  col_stride_mmul_1x2_2x1(&G(0, 0), a_col_stride,
			  &H(0, 0), b_col_stride,
			  &I(0, 0), c_col_stride);

  col_stride_mmul_1x2_2x1(&G(1, 0), a_col_stride,
			  &H(0, 0), b_col_stride,
			  &I(1, 0), c_col_stride);
}

void col_stride_mmul_2x2_2x2(const double *a, int a_col_stride,
			     const double *b, int b_col_stride,
			     double *c, int c_col_stride) {

  col_stride_mmul_2x2_2x1(&G(0, 0), a_col_stride,
			  &H(0, 0), b_col_stride,
			  &I(0, 0), c_col_stride);

  col_stride_mmul_2x2_2x1(&G(0, 0), a_col_stride,
			  &H(0, 1), b_col_stride,
			  &I(0, 1), c_col_stride);
}

#endif
