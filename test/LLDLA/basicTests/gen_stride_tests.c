#include "gen_stride_tests.h"

void general_stride_tests()	{
	int size = 50;
	double a[size];
	double b[size];
	double c[size];
	double c_copy[size];

	reset_values(size, a, b, c, c_copy);

	printf("\n\n---------- START GENERAL STRIDE TESTS ----------\n\n");

	test_buffer_diff(size, c, c_copy, "sanity check");

	// General stride add tests

	simple_add(2, 1, a, 1, 4, c, 1, 4);
	gen_stride_add_2x1(a, 1, 4, c_copy, 1, 4);
	test_buffer_diff(size, c, c_copy, "gen_stride_add_2x1 same row and col stride");

	reset_values(size, a, b, c, c_copy);

	simple_add(2, 1, a, 1, 4, c, 1, 5);
	gen_stride_add_2x1(a, 1, 4, c_copy, 1, 5);
	test_buffer_diff(size, c, c_copy, "gen_stride_add_2x1 same row stride, different col stride");

	reset_values(size, a, b, c, c_copy);

	simple_add(2, 1, a, 1, 4, c, 2, 5);
	gen_stride_add_2x1(a, 1, 4, c_copy, 2, 5);
	test_buffer_diff(size, c, c_copy, "gen_stride_add_2x1 different row and col stride");

	reset_values(size, a, b, c, c_copy);

	simple_add(1, 2, a, 2, 7, c, 2, 7);
	gen_stride_add_1x2(a, 2, 7, c_copy, 2, 7);
	test_buffer_diff(size, c, c_copy, "gen_stride_add_1x2 same row and col stride");

	reset_values(size, a, b, c, c_copy);

	simple_add(1, 2, a, 2, 7, c, 3, 7);
	gen_stride_add_1x2(a, 2, 7, c_copy, 3, 7);
	test_buffer_diff(size, c, c_copy, "gen_stride_add_1x2 different row stride, same col stride");

	reset_values(size, a, b, c, c_copy);

	simple_add(1, 2, a, 2, 7, c, 19, 1);
	gen_stride_add_1x2(a, 2, 7, c_copy, 19, 1);
	test_buffer_diff(size, c, c_copy, "gen_stride_add_1x2 different row and col stride");

	reset_values(size, a, b, c, c_copy);

	simple_add(2, 2, a, 2, 9, c, 2, 9);
	gen_stride_add_2x2(a, 2, 9, c_copy, 2, 9);
	test_buffer_diff(size, c, c_copy, "gen_stride_add_2x2 same row and col stride");

	reset_values(size, a, b, c, c_copy);

	simple_add(2, 2, a, 2, 8, c, 3, 8);
	gen_stride_add_2x2(a, 2, 8, c_copy, 3, 8);
	test_buffer_diff(size, c, c_copy, "gen_stride_add_2x2 different row stride and same col stride");

	reset_values(size, a, b, c, c_copy);

	simple_add(2, 2, a, 2, 6, c, 2, 11);
	gen_stride_add_2x2(a, 2, 6, c_copy, 2, 11);
	test_buffer_diff(size, c, c_copy, "gen_stride_add_2x2 same row stride and different col stride");

	reset_values(size, a, b, c, c_copy);

	simple_add(2, 2, a, 3, 10, c, 5, 16);
	gen_stride_add_2x2(a, 3, 10, c_copy, 5, 16);
	test_buffer_diff(size, c, c_copy, "gen_stride_add_2x2 different row and col stride");

	reset_values(size, a, b, c, c_copy);

	// Multiplication LLDLA tests

	simple_mmul(1, 1, 2, a, 4, 1, b, 4, 1, c, 4, 1);
	gen_stride_mmul_1x2_2x1(a, 4, 1, b, 4, 1, c_copy, 4, 1);
	test_buffer_diff(size, c, c_copy, "gen_stride_mmul_1x2_2x1 same row and col stride");

	reset_values(size, a, b, c, c_copy);

	simple_mmul(1, 1, 2, a, 4, 1, b, 1, 10, c, 7, 9);
	gen_stride_mmul_1x2_2x1(a, 4, 1, b, 1, 10, c_copy, 7, 9);
	test_buffer_diff(size, c, c_copy, "gen_stride_mmul_1x2_2x1 different row and col stride");

	reset_values(size, a, b, c, c_copy);

	simple_mmul(1, 2, 2, a, 5, 2, b, 9, 2, c, 4, 8);
	gen_stride_mmul_1x2_2x2(a, 5, 2, b, 9, 2, c_copy, 4, 8);
	test_buffer_diff(size, c, c_copy, "gen_stride_mmul_1x2_2x2 different row and col stride");

	reset_values(size, a, b, c, c_copy);

	simple_mmul(2, 2, 1, a, 7, 23, b, 6, 1, c, 9, 4);
	gen_stride_mmul_2x1_1x2(a, 7, 23, b, 6, 1, c_copy, 9, 4);
	test_buffer_diff(size, c, c_copy, "gen_stride_mmul_2x1_1x2 different row and col stride");

	reset_values(size, a, b, c, c_copy);

	simple_mmul(2, 1, 2, a, 7, 2, b, 5, 1, c, 19, 14);
	gen_stride_mmul_2x2_2x1(a, 7, 2, b, 5, 1, c_copy, 19, 14);
	test_buffer_diff(size, c, c_copy, "gen_stride_mmul_2x2_2x1 different row and col stride");

	reset_values(size, a, b, c, c_copy);

	simple_mmul(2, 2, 2, a, 11, 3, b, 5, 2, c, 8, 7);
	gen_stride_mmul_2x2_2x2(a, 11, 3, b, 5, 2, c_copy, 8, 7);
	test_buffer_diff(size, c, c_copy, "gen_stride_mmul_2x2_2x2 different row and col stride");

	reset_values(size, a, b, c, c_copy);

	// Scalar multiplication tests

	simple_smul(1, 1, &a[3], c, 5, 2);
	gen_stride_smul_1x1(&a[3], c_copy, 5, 2);
	test_buffer_diff(size, c, c_copy, "gen_stride_smul_1x1");

	reset_values(size, a, b, c, c_copy);

	simple_smul(2, 1, &a[5], c, 11, 4);
	gen_stride_smul_2x1(&a[5], c_copy, 11, 4);
	test_buffer_diff(size, c, c_copy, "gen_stride_smul_2x1");

	reset_values(size, a, b, c, c_copy);

	simple_smul(1, 2, &a[12], c, 3, 9);
	gen_stride_smul_1x2(&a[12], c_copy, 3, 9);
	test_buffer_diff(size, c, c_copy, "gen_stride_smul_1x2");

	reset_values(size, a, b, c, c_copy);

	simple_smul(2, 2, &a[49], c, 15, 1);
	gen_stride_smul_2x2(&a[49], c_copy, 15, 1);
	test_buffer_diff(size, c, c_copy, "gen_stride_smul_2x2");

	printf("\n\n---------- END GENERAL STRIDE TESTS ----------\n\n");

	return;
}