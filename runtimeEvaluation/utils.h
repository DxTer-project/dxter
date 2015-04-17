#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <pmmintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

void* alloc(size_t size);

void* alloc_aligned_16(size_t size);
void* alloc_aligned_32(size_t size);

void pack_double(double* a, double* b,
		 int a_rows, int a_cols,
		 int a_row_stride, int a_col_stride,
		 int b_rows, int b_cols,
		 int b_row_stride, int b_col_stride);

void pack_float(float* a, float* b,
		int a_rows, int a_cols,
		int a_row_stride, int a_col_stride,
		int b_rows, int b_cols,
		int b_row_stride, int b_col_stride);

void unpack_double(double* a, double* b,
		   int a_row_stride, int a_col_stride,
		   int b_rows, int b_cols,
		   int b_row_stride, int b_col_stride);

void unpack_float(float* a, float* b,
		  int a_row_stride, int a_col_stride,
		  int b_rows, int b_cols,
		  int b_row_stride, int b_col_stride);

void set_to_zero_double(double* a,
			int m, int n,
			int a_row_stride, int a_col_stride);

void set_to_zero_float(float* a,
		       int m, int n,
		       int a_row_stride, int a_col_stride);

void rand_doubles(int size, double* rands);

void rand_floats(int size, float *rands);

void simple_add(int m, int n,
	double* a, int a_row_stride, int a_col_stride,
	double* b, int b_row_stride, int b_col_stride);

void simple_mmul(int m, int n, int k,
	double* a, int a_row_stride, int a_col_stride,
	double* b, int b_row_stride, int b_col_stride,
	double* c, int c_row_stride, int c_col_stride);

void simple_smul(int m, int n,
	double* scalar,
	double* a, int a_row_stride, int a_col_stride);

void simple_add_float(int m, int n,
	float *a, int a_row_stride, int a_col_stride,
	float *b, int b_row_stride, int b_col_stride);

void simple_mmul_float(int m, int n, int k,
	float *a, int a_row_stride, int a_col_stride,
	float *b, int b_row_stride, int b_col_stride,
	float *c, int c_row_stride, int c_col_stride);

void simple_smul_float(int m, int n,
	float *scalar,
	float *a, int a_row_stride, int a_col_stride);

void copy_double(int m, int n,
		 double* a,
		 int a_row_stride, int a_col_stride,
		 double* b,
		 int b_row_stride, int b_col_stride);

void copy_float(int m, int n,
		 float* a,
		 int a_row_stride, int a_col_stride,
		 float* b,
		 int b_row_stride, int b_col_stride);

double diff_buffer(int size, double* buf1, double* buf2);

float diff_buffer_float(int size, float *buf1, float *buf2);

void copy_buffer(int size, double* src, double* dest);

void copy_buffer_float(int size, float *src, float *dest);

void print_buffer(int size, double* buf);

void test_buffer_diff(int size, double* a, double* b, char *test_name);

void test_buffer_diff_float(int size, float *a, float *b, char *test_name);

void test_mats_diff(int m, int n, double* a, int a_row_stride, int a_col_stride, double* b, int b_row_stride, int b_col_stride, char *test_name);

void print_mat(int m, int n, double* a, int a_row_stride, int a_col_stride);
void print_mat_float(int m, int n, float* a, int a_row_stride, int a_col_stride);

void reset_values(int size, double* a, double* b, double* c, double* c_copy);

unsigned long long rdtsc();

#endif
