#ifndef GEN_STRIDE_LLDLA_PRIMITIVES_H
#define GEN_STRIDE_LLDLA_PRIMITIVES_H

// Add operations

inline void gen_stride_add_2x1(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride);

inline void gen_stride_add_1x2(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride);

inline void gen_stride_add_2x2(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride);

// Matrix matrix multiplication operations

inline void gen_stride_mmul_1x2_2x1(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride);

inline void gen_stride_mmul_1x2_2x2(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride);

inline void gen_stride_mmul_2x1_1x2(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride);

inline void gen_stride_mmul_2x2_2x1(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride);

inline void gen_stride_mmul_2x2_2x2(
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride);

// Scalar multiply operations

inline void gen_stride_smul_1x1(
	double *scalar,
	double *a, int a_row_stride, int a_col_stride);

inline void gen_stride_smul_2x1(
	double *scalar,
	double *a, int a_row_stride, int a_col_stride);

inline void gen_stride_smul_1x2(
	double *scalar,
	double *a, int a_row_stride, int a_col_stride);

inline void gen_stride_smul_2x2(
	double *scalar,
	double *a, int a_row_stride, int a_col_stride);

#endif
