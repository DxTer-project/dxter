#ifndef COL_STRIDE_LLDLA_PRIMITIVES_H
#define COL_STRIDE_LLDLA_PRIMITIVES_H

#include "utils.h"

#define A(i, j) a[ (i) + (j)*a_col_stride ]
#define B(i, j) b[ (i) + (j)*b_col_stride ]
#define C(i, j) c[ (i) + (j)*c_col_stride ]

// Addition kernels

void col_stride_add_2x1(double *a, int a_col_stride,
			double *b, int b_col_stride) {
  v2df_t b_00_10;
  b_00_10.v = VEC_PD_LOAD(&B(0, 0));
  b_00_10.v += VEC_PD_LOAD(&A(0, 0));
  VEC_PD_STORE(&B(0, 0), b_00_10.v);
}

void col_stride_add_1x2(double *a, int a_col_stride,
			double *b, int b_col_stride) {
  v2df_t a_00_01, b_00_01;
  a_00_01.v = VEC_D_LOAD(&A(0, 0), &A(0, 1));
  b_00_01.v = VEC_D_LOAD(&B(0, 0), &B(0, 1));
  b_00_01.v += a_00_01.v;
  B(0, 0) = b_00_01.d[0];
  B(0, 1) = b_00_01.d[1];
}

void col_stride_add_2x2(double *a, int a_col_stride,
			double *b, int b_col_stride) {
  col_stride_add_2x1(&A(0, 0), a_col_stride, &B(0, 0), b_col_stride);
  col_stride_add_2x1(&A(0, 1), a_col_stride, &B(0, 1), b_col_stride);
}

// Scalar multiplication kernels

void col_stride_smul_2x1(const double *scalar,
			 double *a, int a_col_stride) {

  v2df_t scal, a_00_10;
  
  scal.v = VEC_DUP_LOAD(scalar);
  a_00_10.v = VEC_PD_LOAD(&A(0, 0));
  a_00_10.v *= scal.v;
  VEC_PD_STORE(&A(0, 0), a_00_10.v);
}

void col_stride_smul_1x2(const double *scalar,
			 double *a, int a_col_stride) {
  v2df_t scal, a_00_01;
  scal.v = VEC_DUP_LOAD(scalar);
  a_00_01.v = VEC_D_LOAD(&A(0, 0), &A(0, 1));
  a_00_01.v = a_00_01.v * scal.v;
  A(0, 0) = a_00_01.d[0];
  A(0, 1) = a_00_01.d[1];
}

void col_stride_smul_2x2(const double *scalar,
			 double *a, int a_col_stride) {
  col_stride_smul_2x1(scalar, &A(0, 0), a_col_stride);
  col_stride_smul_2x1(scalar, &A(0, 1), a_col_stride);
}

// Matrix multiply primitives

void col_stride_mmul_1x2_2x1(double *a, int a_col_stride,
			     double *b, int b_col_stride,
			     double *c, int c_col_stride) {
  v2df_t a_00_01, b_00_10;
  
  a_00_01.v = VEC_D_LOAD(&A(0, 0), &A(0, 1));
  b_00_10.v = VEC_PD_LOAD(&B(0, 0));
  a_00_01.v *= b_00_10.v;
  C(0, 0) += a_00_01.d[0];
  C(0, 0) += a_00_01.d[1];
}

void col_stride_mmul_2x1_1x2(double *a, int a_col_stride,
			     double *b, int b_col_stride,
			     double *c, int c_col_stride) {

  v2df_t a_00_10, b_00, b_01;

  a_00_10.v = VEC_PD_LOAD(&A(0, 0));
  b_00.v = VEC_DUP_LOAD(&B(0, 0));
  b_01.v = VEC_DUP_LOAD(&B(0, 1));

  b_00.v *= a_00_10.v;
  b_01.v *= a_00_10.v;

  C(0, 0) += b_00.d[0];
  C(1, 0) += b_00.d[1];
  C(0, 1) += b_01.d[0];
  C(1, 1) += b_01.d[1];
}

void col_stride_mmul_1x2_2x2(double *a, int a_col_stride,
			     double *b, int b_col_stride,
			     double *c, int c_col_stride) {

  col_stride_mmul_1x2_2x1(&A(0, 0), a_col_stride,
			  &B(0, 0), b_col_stride,
			  &C(0, 0), c_col_stride);

  col_stride_mmul_1x2_2x1(&A(0, 0), a_col_stride,
			  &B(0, 1), b_col_stride,
			  &C(0, 1), c_col_stride);
}

void col_stride_mmul_2x2_2x1(double *a, int a_col_stride,
			     double *b, int b_col_stride,
			     double *c, int c_col_stride) {

  col_stride_mmul_1x2_2x1(&A(0, 0), a_col_stride,
			  &B(0, 0), b_col_stride,
			  &C(0, 0), c_col_stride);

  col_stride_mmul_1x2_2x1(&A(1, 0), a_col_stride,
			  &B(0, 0), b_col_stride,
			  &C(1, 0), c_col_stride);
}

void col_stride_mmul_2x2_2x2(double *a, int a_col_stride,
			     double *b, int b_col_stride,
			     double *c, int c_col_stride) {

  col_stride_mmul_2x2_2x1(&A(0, 0), a_col_stride,
			  &B(0, 0), b_col_stride,
			  &C(0, 0), c_col_stride);

  col_stride_mmul_2x2_2x1(&A(0, 0), a_col_stride,
			  &B(0, 1), b_col_stride,
			  &C(0, 1), c_col_stride);
}

#endif
