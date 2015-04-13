#include "utils.h"

#define A(i,j) a[ (i)*a_row_stride + (j)*a_col_stride ]
#define B(i,j) b[ (i)*b_row_stride + (j)*b_col_stride ]
#define C(i,j) c[ (i)*c_row_stride + (j)*c_col_stride ]

void *alloc(size_t size)	{
	return (void *) malloc(size);
}

void pack_double(double* a, double* b,
		 int a_rows, int a_cols,
		 int a_row_stride, int a_col_stride,
		 int b_rows, int b_cols,
		 int b_row_stride, int b_col_stride) {
  int i, j;
  for (i = 0; i < b_rows; i++) {
    for (j = 0; j < b_cols; j++) {
      if (i < a_rows && j < a_cols) {
	B(i, j) = A(i, j);
      } else {
	B(i, j) = 0;
      }
    }
  }
}

void pack_float(float* a, float* b,
		int a_rows, int a_cols,
		int a_row_stride, int a_col_stride,
		int b_rows, int b_cols,
		int b_row_stride, int b_col_stride) {
  int i, j;
  for (i = 0; i < b_rows; i++) {
    for (j = 0; j < b_cols; j++) {
      if (i < a_rows && j < a_cols) {
	B(i, j) = A(i, j);
      } else {
	B(i, j) = 0;
      }
    }
  }
}

void unpack_double(double* a, double* b,
		   int a_row_stride, int a_col_stride,
		   int b_rows, int b_cols,
		   int b_row_stride, int b_col_stride) {
  int i, j;
  for (i = 0; i < b_rows; i++) {
    for (j = 0; j < b_cols; j++) {
      B(i, j) = A(i, j);
    }
  }
}

void unpack_float(float* a, float* b,
		  int a_row_stride, int a_col_stride,
		  int b_rows, int b_cols,
		  int b_row_stride, int b_col_stride) {
  int i, j;
  for (i = 0; i < b_rows; i++) {
    for (j = 0; j < b_cols; j++) {
      B(i, j) = A(i, j);
    }
  }
}

void set_to_zero_double(double* a,
			int m, int n,
			int a_row_stride, int a_col_stride) {
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      A(i, j) = 0;
    }
  }
}

void set_to_zero_float(float* a,
		       int m, int n,
		       int a_row_stride, int a_col_stride) {
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      A(i, j) = 0;
    }
  }
}

void *alloc_aligned_16(size_t size) {
  void *memptr;
  // Horrible error handling
  if (posix_memalign(&memptr, 16, size))	{
  	printf("FAILED TO ALLOCATE\n");
  	return (void *) -1;
  }
  return (void *) memptr;
}

void *alloc_aligned_32(size_t size) {
  void *memptr;
  // Horrible error handling
  if (posix_memalign(&memptr, 32, size))	{
  	printf("FAILED TO ALLOCATE\n");
  	return (void *) -1;
  }
  return (void *) memptr;
}

void copy_double(int m, int n,
		 double* a,
		 int a_row_stride, int a_col_stride,
		 double* b,
		 int b_row_stride, int b_col_stride) {
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      B(i, j) = A(i, j);
    }
  }
}

void copy_float(int m, int n,
		 float* a,
		 int a_row_stride, int a_col_stride,
		 float* b,
		 int b_row_stride, int b_col_stride) {
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      B(i, j) = A(i, j);
    }
  }
}

void rand_doubles(int size, double *rands) {
  int i;
  for (i = 0; i < size; i++) {
    rands[i] = rand() % 10;
  }
}

void rand_floats(int size, float *rands) {
  int i;
  for (i = 0; i < size; i++) {
    rands[i] = (float) (rand() % 10);
  }
}

void simple_add(int m, int n,
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride) {

  int i, j;
  for (i = 0; i < m; i++)	{
  	for (j = 0; j < n; j++)	{
  		B(i, j) += A(i, j);
  	}
  }
}

void simple_mmul(int m, int n, int k,
	double *a, int a_row_stride, int a_col_stride,
	double *b, int b_row_stride, int b_col_stride,
	double *c, int c_row_stride, int c_col_stride)	{

	int i, j, p;
	for (i = 0; i < m; i++)	{
		for (j = 0; j < n; j++)	{
			for (p = 0; p < k; p++)	{
				C(i, j) += A(i, p) * B(p, j);
			}
		}
	}
}

void simple_smul(int m, int n,
	double *scalar,
	double *a, int a_row_stride, int a_col_stride)	{

	int i, j;
	for (i = 0; i < m; i++)	{
		for (j = 0; j < n; j++)	{
			A(i, j) = A(i, j) * *scalar;
		}
	}
}

void simple_add_float(int m, int n,
		      float *a, int a_row_stride, int a_col_stride,
		      float *b, int b_row_stride, int b_col_stride) {
  int i, j;
  for (i = 0; i < m; i++)	{
    for (j = 0; j < n; j++)	{
      B(i, j) += A(i, j);
    }
  }
  return;
}

void simple_mmul_float(int m, int n, int k,
		       float *a, int a_row_stride, int a_col_stride,
		       float *b, int b_row_stride, int b_col_stride,
		       float *c, int c_row_stride, int c_col_stride) {

  int i, j, p;
  for (i = 0; i < m; i++)	{
    for (j = 0; j < n; j++)	{
      for (p = 0; p < k; p++)	{
	C(i, j) += A(i, p) * B(p, j);
      }
    }
  }
  return;
}

void simple_smul_float(int m, int n,
		       float *scalar,
		       float *a, int a_row_stride, int a_col_stride) {
  int i, j;
  for (i = 0; i < m; i++)	{
    for (j = 0; j < n; j++)	{
      A(i, j) = A(i, j) * *scalar;
    }
  }
  return;
}

double diff_buffer(int size, double *buf1, double *buf2) {
  int i;
  double diff = 0.0;
  for (i = 0; i < size; i++) {
    diff += fabs(buf1[i] - buf2[i]);
  }
  return diff;
}

float diff_buffer_float(int size, float *buf1, float *buf2) {
  int i;
  float diff = 0.0;
  for (i = 0; i < size; i++) {
    diff += fabs(buf1[i] - buf2[i]);
  }
  return diff;
}

double diff_mats(int m, int n, double *a, int a_row_stride, int a_col_stride, double *b, int b_row_stride, int b_col_stride)	{
	int i, j;
	double diff = 0.0;
	for (i = 0; i < m; i++)	{
		for (j = 0; j < n; j++)	{
			diff += fabs(A(i, j) - B(i, j));
		}
	}
	return diff;
}

void test_mats_diff(int m, int n, double *a, int a_row_stride, int a_col_stride, double *b, int b_row_stride, int b_col_stride, char *test_name)	{
	double diff = diff_mats(m, n, a, a_row_stride, a_col_stride, b, b_row_stride, b_col_stride);
	// Should really have a tolerance but for now this will do
	if (diff != 0.0)	{
		printf("\n\nERROR in %s: diff = %f\n\n", test_name, diff);
	} else {
		printf("%s PASSED\n", test_name);
	}
}

void copy_buffer(int size, double *src, double *dest) {
  int i;
  for (i = 0; i < size; i++) {
    dest[i] = src[i];
  }
}

void copy_buffer_float(int size, float *src, float *dest) {
  int i;
  for (i = 0; i < size; i++) {
    dest[i] = src[i];
  }
}

void print_buffer(int size, double *buf)	{
	int i;
	for (i = 0; i < size; i++)	{
		printf("%f ", buf[i]);
	}
	printf("\n");
}

void test_buffer_diff(int size, double *a, double *b, char *test_name)	{
	double diff = diff_buffer(size, a, b);
	// Should really have a tolerance but for now this will do
	if (diff != 0.0)	{
		printf("\n\nERROR in %s: diff = %f\n\n", test_name, diff);
	} else {
		printf("%s PASSED\n", test_name);
	}
}

void test_buffer_diff_float(int size, float *a, float *b, char *test_name)	{
	float diff = diff_buffer_float(size, a, b);
	// Should really have a tolerance but for now this will do
	if (diff != 0.0)	{
		printf("\n\nERROR in %s: diff = %f\n\n", test_name, diff);
	} else {
		printf("%s PASSED\n", test_name);
	}
}

void print_mat(int m, int n, double *a, int a_row_stride, int a_col_stride)	{
	int i, j;
	for (i = 0; i < m; i++)	{
		for (j = 0; j < n; j++)	{
			printf("%f ", A(i, j));
		}
		printf("\n");
	}
	printf("\n");
}

void reset_values(int size, double *a, double *b, double *c, double *c_copy)	{	rand_doubles(size, a);	rand_doubles(size, b);	rand_doubles(size, c);	copy_buffer(size, c, c_copy);}
