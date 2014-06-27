#ifndef UTILS_H
#define UTILS_H

//#include <malloc.h>
#include <math.h>
#include <pmmintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define MUVALUE 2
#define VEC_PD_LOAD(ptr) _mm_load_pd((ptr))									// Load array of 2 doubles to vector register
#define VEC_DUP_LOAD(ptr) _mm_loaddup_pd((ptr))								// Load one double into upper and lower halves of vector register
#define VEC_D_LOAD(ptr1, ptr2) _mm_loadh_pd(_mm_load_sd((ptr1)), (ptr2));	// Load 2 doubles at different memory locations into the 2 halves of a vector register

#define VEC_PD_STORE(ptr, vec_reg) _mm_store_pd((ptr), (vec_reg))			// Store vector register contents in array of 2 doubles

typedef union	{
	__m128d v;
	double d[2];
} v2df_t;

void *alloc(size_t size);

void *alloc_aligned_16(size_t size);

void rand_doubles(int size, double *rands);

void simple_add(int m, int n,
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride);

void simple_mmul(int m, int n, int k,
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride);

void simple_smul(int m, int n,
	double *scalar,
	double *a, int a_row_stride, int a_col_stride);

double diff_buffer(int size, double *buf1, double *buf2);

void copy_buffer(int size, double *src, double *dest);

void print_buffer(int size, double *buf);

void test_buffer_diff(int size, double *a, double *b, char *test_name);

void test_mats_diff(int m, int n, double *a, int a_row_stride, int a_col_stride, double *b, int b_row_stride, int b_col_stride, char *test_name);

void print_mat(int m, int n, double *a, int a_row_stride, int a_col_stride);

void reset_values(int size, double *a, double *b, double *c, double *c_copy);

#endif
