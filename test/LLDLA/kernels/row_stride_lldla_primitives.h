#ifndef ROW_STRIDE_LLDLA_PRIMITIVES_H
#define ROW_STRIDE_LLDLA_PRIMITIVES_H

#include "utils.h"

#define A(i,j) a[ (i)*a_row_stride + j ]
#define B(i,j) b[ (i)*b_row_stride + j ]
#define C(i,j) c[ (i)*c_row_stride + j ]

// Add operations

inline void row_stride_add_2x1(
	double *a, int a_row_stride,
	double *b, int b_row_stride)	{

	v2df_t
		a_00_10_vreg,
		b_00_10_vreg;

	a_00_10_vreg.v = VEC_D_LOAD( &A(0, 0), &A(1, 0) );
	b_00_10_vreg.v = VEC_D_LOAD( &B(0, 0), &B(1, 0) );
	b_00_10_vreg.v += a_00_10_vreg.v;

	B(0, 0) = b_00_10_vreg.d[0];
	B(1, 0) = b_00_10_vreg.d[1];
}

inline void row_stride_add_1x2(
	double *a, int a_row_stride,
	double *b, int b_row_stride)	{

	v2df_t
		b_00_01_vreg;

	b_00_01_vreg.v = VEC_PD_LOAD( &B(0, 0) );
	b_00_01_vreg.v += VEC_PD_LOAD( &A(0, 0) );
	VEC_PD_STORE( &B(0, 0), b_00_01_vreg.v );
}

inline void row_stride_add_2x2(
	double *a, int a_row_stride,
	double *b, int b_row_stride)	{

	row_stride_add_1x2(
		&A(0, 0), a_row_stride,
		&B(0, 0), b_row_stride);

	row_stride_add_1x2(
		&A(1, 0), a_row_stride,
		&B(1, 0), b_row_stride);
}

// Matrix multiply operations

inline void row_stride_mmul_1x2_2x1(
	double *a, int a_row_stride,
	double *b, int b_row_stride,
	double *c, int c_row_stride)	{

	v2df_t
		a_00_01_vreg,
		b_00_10_vreg;

	a_00_01_vreg.v = VEC_PD_LOAD(&A(0, 0));
	b_00_10_vreg.v = VEC_D_LOAD(&B(0, 0), &B(1, 0));
	b_00_10_vreg.v *= a_00_01_vreg.v;

	C(0, 0) += b_00_10_vreg.d[0];
	C(0, 0) += b_00_10_vreg.d[1];
}

inline void row_stride_mmul_1x2_2x2(
	double *a, int a_row_stride,
	double *b, int b_row_stride,
	double *c, int c_row_stride)	{

	v2df_t
		a_00_vreg, a_01_vreg,
		b_00_01_vreg, b_10_11_vreg,
		c_00_01_vreg;

	b_00_01_vreg.v = VEC_PD_LOAD(&B(0, 0));
	b_10_11_vreg.v = VEC_PD_LOAD(&B(1, 0));

	a_00_vreg.v = VEC_DUP_LOAD(&A(0, 0));
	a_01_vreg.v = VEC_DUP_LOAD(&A(0, 1));

	b_00_01_vreg.v *= a_00_vreg.v;
	b_10_11_vreg.v *= a_01_vreg.v;

	c_00_01_vreg.v = VEC_PD_LOAD(&C(0, 0));
	c_00_01_vreg.v += b_00_01_vreg.v + b_10_11_vreg.v;

	VEC_PD_STORE(&C(0, 0), c_00_01_vreg.v);
}

inline void row_stride_mmul_2x1_1x2(
	double *a, int a_row_stride,
	double *b, int b_row_stride,
	double *c, int c_row_stride)	{

	v2df_t
		a_00_vreg, a_10_vreg,
		b_00_01_vreg,
		c_00_01_vreg, c_10_11_vreg;

	b_00_01_vreg.v = VEC_PD_LOAD(&B(0, 0));

	a_00_vreg.v = VEC_DUP_LOAD(&A(0, 0));
	a_10_vreg.v = VEC_DUP_LOAD(&A(1, 0));

	c_00_01_vreg.v = VEC_PD_LOAD(&C(0, 0));
	c_10_11_vreg.v = VEC_PD_LOAD(&C(1, 0));

	c_00_01_vreg.v += a_00_vreg.v * b_00_01_vreg.v;
	c_10_11_vreg.v += a_10_vreg.v * b_00_01_vreg.v;

	VEC_PD_STORE(&C(0, 0), c_00_01_vreg.v);
	VEC_PD_STORE(&C(1, 0), c_10_11_vreg.v);
}

inline void row_stride_mmul_2x2_2x1(
	double *a, int a_row_stride,
	double *b, int b_row_stride,
	double *c, int c_row_stride)	{

	v2df_t
		a_00_01_vreg, a_10_11_vreg,
		b_00_10_vreg,
		c_00_10_vreg;

	a_00_01_vreg.v = VEC_PD_LOAD(&A(0, 0));
	a_10_11_vreg.v = VEC_PD_LOAD(&A(1, 0));

	b_00_10_vreg.v = VEC_D_LOAD(&B(0, 0), &B(1, 0));

	a_00_01_vreg.v *= b_00_10_vreg.v;
	a_10_11_vreg.v *= b_00_10_vreg.v;

	C(0, 0) += a_00_01_vreg.d[0] + a_00_01_vreg.d[1];
	C(1, 0) += a_10_11_vreg.d[0] + a_10_11_vreg.d[1];
}

inline void row_stride_mmul_2x2_2x2(
	double *a, int a_row_stride,
	double *b, int b_row_stride,
	double *c, int c_row_stride)	{

	v2df_t
		a_00_vreg, a_01_vreg, a_10_vreg, a_11_vreg,
		b_00_10_vreg, b_10_11_vreg,
		c_00_01_vreg, c_10_11_vreg;

	a_00_vreg.v = VEC_DUP_LOAD(&A(0, 0));
	a_10_vreg.v = VEC_DUP_LOAD(&A(1, 0));
	a_01_vreg.v = VEC_DUP_LOAD(&A(0, 1));
	a_11_vreg.v = VEC_DUP_LOAD(&A(1, 1));

	b_00_10_vreg.v = VEC_PD_LOAD(&B(0, 0));
	b_10_11_vreg.v = VEC_PD_LOAD(&B(1, 0));

	c_00_01_vreg.v = VEC_PD_LOAD(&C(0, 0));
	c_00_01_vreg.v += a_00_vreg.v * b_00_10_vreg.v + a_01_vreg.v * b_10_11_vreg.v;
	c_10_11_vreg.v = VEC_PD_LOAD(&C(1, 0));
	c_10_11_vreg.v += a_10_vreg.v * b_00_10_vreg.v + a_11_vreg.v * b_10_11_vreg.v;

	VEC_PD_STORE(&C(0, 0), c_00_01_vreg.v);
	VEC_PD_STORE(&C(1, 0), c_10_11_vreg.v);
}

// Scalar multiply operations

inline void row_stride_smul_1x1(
	double *scalar,
	double *a, int a_row_stride)	{

	A(0, 0) =  A(0, 0) * *scalar;
}

inline void row_stride_smul_2x1(
	double *scalar,
	double *a, int a_row_stride)	{

	v2df_t
		scalar_vreg,
		a_00_10_vreg;

	scalar_vreg.v = VEC_DUP_LOAD(scalar);

	a_00_10_vreg.v = VEC_D_LOAD(&A(0, 0), &A(1, 0));
	a_00_10_vreg.v *= scalar_vreg.v;

	A(0, 0) = a_00_10_vreg.d[0];
	A(1, 0) = a_00_10_vreg.d[1];
}

inline void row_stride_smul_1x2(
	double *scalar,
	double *a, int a_row_stride)	{

	v2df_t
		scalar_vreg,
		a_00_10_vreg;

	scalar_vreg.v = VEC_DUP_LOAD(scalar);

	a_00_10_vreg.v = VEC_PD_LOAD(&A(0, 0));
	a_00_10_vreg.v *= scalar_vreg.v;

	VEC_PD_STORE(&A(0, 0), a_00_10_vreg.v);
}

inline void row_stride_smul_2x2(
	double *scalar,
	double *a, int a_row_stride)	{

	v2df_t
		scalar_vreg,
		a_00_10_vreg, a_10_11_vreg;

	scalar_vreg.v = VEC_DUP_LOAD(scalar);

	a_00_10_vreg.v = VEC_PD_LOAD(&A(0, 0));
	a_10_11_vreg.v = VEC_PD_LOAD(&A(1, 0));
	a_00_10_vreg.v *= scalar_vreg.v;
	a_10_11_vreg.v *= scalar_vreg.v;

	VEC_PD_STORE(&A(0, 0), a_00_10_vreg.v);
	VEC_PD_STORE(&A(1, 0), a_10_11_vreg.v);
}
#endif
