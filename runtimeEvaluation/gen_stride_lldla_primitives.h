#ifndef GEN_STRIDE_LLDLA_PRIMITIVES_H
#define GEN_STRIDE_LLDLA_PRIMITIVES_H

#include "utils.h"

#define D(i,j) a[ (i)*a_row_stride + (j)*a_col_stride ]
#define E(i,j) b[ (i)*b_row_stride + (j)*b_col_stride ]
#define F(i,j) c[ (i)*c_row_stride + (j)*c_col_stride ]

// Addition operations

inline void gen_stride_add_2x1(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride)	{
	
	v2df_t
		a_00_10_vreg,
		b_00_10_vreg;
	
	a_00_10_vreg.v = VEC_D_LOAD( &D(0, 0), &D(1, 0) );
	b_00_10_vreg.v = VEC_D_LOAD( &E(0, 0), &E(1, 0) );
	b_00_10_vreg.v += a_00_10_vreg.v;

	E(0, 0) = b_00_10_vreg.d[0];
	E(1, 0) = b_00_10_vreg.d[1];
}

inline void gen_stride_add_1x2(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride)	{

	v2df_t
		a_00_01_vreg,
		b_00_01_vreg;
	
	a_00_01_vreg.v = VEC_D_LOAD(&D(0, 0), &D(0, 1));
	b_00_01_vreg.v = VEC_D_LOAD(&E(0, 0), &E(0, 1));
	b_00_01_vreg.v += a_00_01_vreg.v;

	E(0, 0) = b_00_01_vreg.d[0];
	E(0, 1) = b_00_01_vreg.d[1];
}

inline void gen_stride_add_2x2(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride)	{

	gen_stride_add_2x1(&D(0, 0), a_row_stride, a_col_stride, &E(0, 0), b_row_stride, b_col_stride);
	gen_stride_add_2x1(&D(0, 1), a_row_stride, a_col_stride, &E(0, 1), b_row_stride, b_col_stride);
}

// Matrix multiply operations

inline void gen_stride_mmul_1x2_2x1(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride)	{

	v2df_t
		a_00_01_vreg,
		b_00_10_vreg;

	a_00_01_vreg.v = VEC_D_LOAD(&D(0, 0), &D(0, 1));
	b_00_10_vreg.v = VEC_D_LOAD(&E(0, 0), &E(1, 0));
	b_00_10_vreg.v *= a_00_01_vreg.v;

	F(0, 0) += b_00_10_vreg.d[0];
	F(0, 0) += b_00_10_vreg.d[1];
}

inline void gen_stride_mmul_1x2_2x2(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride)	{

	gen_stride_mmul_1x2_2x1(
		&D(0, 0), a_row_stride, a_col_stride,
		&E(0, 0), b_row_stride, b_col_stride,
		&F(0, 0), c_row_stride, c_col_stride);

	gen_stride_mmul_1x2_2x1(
		&D(0, 0), a_row_stride, a_col_stride,
		&E(0, 1), b_row_stride, b_col_stride,
		&F(0, 1), c_row_stride, c_col_stride);
}

inline void gen_stride_mmul_2x1_1x2(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride)	{

	v2df_t
		a_00_10_vreg,
		b_00_vreg, b_01_vreg,
		c_00_10_vreg, c_01_11_vreg;

	a_00_10_vreg.v = VEC_D_LOAD(&D(0, 0), &D(1, 0));
	
	b_00_vreg.v = VEC_DUP_LOAD(&E(0, 0));
	b_01_vreg.v = VEC_DUP_LOAD(&E(0, 1));

	c_00_10_vreg.v = a_00_10_vreg.v * b_00_vreg.v;
	c_01_11_vreg.v = a_00_10_vreg.v * b_01_vreg.v;

	F(0, 0) += c_00_10_vreg.d[0];
	F(1, 0) += c_00_10_vreg.d[1];
	F(0, 1) += c_01_11_vreg.d[0];
	F(1, 1) += c_01_11_vreg.d[1];
}

inline void gen_stride_mmul_2x2_2x1(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride)	{

	gen_stride_mmul_1x2_2x1(
		&D(0, 0), a_row_stride, a_col_stride,
		&E(0, 0), b_row_stride, b_col_stride,
		&F(0, 0), c_row_stride, c_col_stride);

	gen_stride_mmul_1x2_2x1(
		&D(1, 0), a_row_stride, a_col_stride,
		&E(0, 0), b_row_stride, b_col_stride,
		&F(1, 0), c_row_stride, c_col_stride);
}

inline void gen_stride_mmul_2x2_2x2(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride)	{

	gen_stride_mmul_2x2_2x1(
		&D(0, 0), a_row_stride, a_col_stride,
		&E(0, 0), b_row_stride, b_col_stride,
		&F(0, 0), c_row_stride, c_col_stride);

	gen_stride_mmul_2x2_2x1(
		&D(0, 0), a_row_stride, a_col_stride,
		&E(0, 1), b_row_stride, b_col_stride,
		&F(0, 1), c_row_stride, c_col_stride);
}

// Scalar multiplication operations

inline void gen_stride_smul_1x1(
	double *scalar,
	double *a, int a_row_stride, int a_col_stride)	{

	D(0, 0) =  D(0, 0) * *scalar;
}

inline void gen_stride_smul_2x1(
	double *scalar,
	double *a, int a_row_stride, int a_col_stride)	{

	v2df_t
		a_00_10_vreg,
		scalar_vreg;

	scalar_vreg.v = VEC_DUP_LOAD(scalar);

	a_00_10_vreg.v = VEC_D_LOAD(&D(0, 0), &D(1, 0));
	a_00_10_vreg.v *= scalar_vreg.v;

	D(0, 0) = a_00_10_vreg.d[0];
	D(1, 0) = a_00_10_vreg.d[1];
}

inline void gen_stride_smul_1x2(
	double *scalar,
	double *a, int a_row_stride, int a_col_stride)	{

	v2df_t
		a_00_01_vreg,
		scalar_vreg;

	scalar_vreg.v = VEC_DUP_LOAD(scalar);

	a_00_01_vreg.v = VEC_D_LOAD(&D(0, 0), &D(0, 1));
	a_00_01_vreg.v *= scalar_vreg.v;

	D(0, 0) = a_00_01_vreg.d[0];
	D(0, 1) = a_00_01_vreg.d[1];
}

inline void gen_stride_smul_2x2(
	double *scalar,
	double *a, int a_row_stride, int a_col_stride)	{

	gen_stride_smul_2x1(
		scalar,
		&D(0, 0), a_row_stride, a_col_stride);

	gen_stride_smul_2x1(
		scalar,
		&D(0, 1), a_row_stride, a_col_stride);
}

#endif
