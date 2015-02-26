#include "row_stride_lldla_primitives.h"
#include "utils.h"

void row_stride_tests()	{

	int size = 50;
	double *a = alloc_aligned_16(size * sizeof(double));
	double *b = alloc_aligned_16(size * sizeof(double));
	double *c = alloc_aligned_16(size * sizeof(double));
	double *c_copy = alloc_aligned_16(size * sizeof(double));

	reset_values(size, a, b, c, c_copy);

	printf("\n\n---------- START ROW STRIDE TESTS ----------\n\n");

	test_buffer_diff(size, c, c_copy, "sanity check");

	// Row stride add tests

	simple_add(2, 1, a, 4, 1, c, 7, 1);
	row_stride_add_2x1(a, 4, c_copy, 7);
	test_buffer_diff(size, c, c_copy, "row_stride_add_2x1");

	reset_values(size, a, b, c, c_copy);

	simple_add(1, 2, a, 10, 1, c, 4, 1);
	row_stride_add_1x2(a, 10, c_copy, 4);
	test_buffer_diff(size, c, c_copy, "row_stride_add_1x2");

	reset_values(size, a, b, c, c_copy);

	simple_add(2, 2, a, 4, 1, c, 2, 1);
	row_stride_add_2x2(a, 4, c_copy, 2);
	test_buffer_diff(size, c, c_copy, "row_stride_add_2x2");

	// Row stride matrix multiplication tests

	reset_values(size, a, b, c, c_copy);

	simple_mmul(1, 1, 2, a, 4, 1, b, 2, 1, c, 6, 1);
	row_stride_mmul_1x2_2x1(a, 4, b, 2, c_copy, 6);
	test_buffer_diff(size, c, c_copy, "row_stride_mmul_1x2_2x1");

	reset_values(size, a, b, c, c_copy);

	simple_mmul(2, 2, 1, a, 4, 1, b, 4, 1, c, 6, 1);
	row_stride_mmul_2x1_1x2(a, 4, b, 4, c_copy, 6);
	test_buffer_diff(size, c, c_copy, "row_stride_mmul_2x1_1x2");

	reset_values(size, a, b, c, c_copy);

	simple_mmul(1, 1, 2, a, 4, 1, b, 4, 1, c, 6, 1);
	row_stride_mmul_1x2_2x1(a, 4, b, 4, c_copy, 6);
	test_buffer_diff(size, c, c_copy, "row_stride_mmul_1x2_2x1");

	reset_values(size, a, b, c, c_copy);

	simple_mmul(2, 1, 2, a, 6, 1, b, 4, 1, c, 8, 1);
	row_stride_mmul_2x2_2x1(a, 6, b, 4, c_copy, 8);
	test_buffer_diff(size, c, c_copy, "row_stride_mmul_2x2_2x1");

	reset_values(size, a, b, c, c_copy);

	simple_mmul(2, 2, 2, a, 8, 1, b, 2, 1, c, 8, 1);
	row_stride_mmul_2x2_2x2(a, 8, b, 2, c_copy, 8);
	test_buffer_diff(size, c, c_copy, "row_stride_mmul_2x2_2x2");

	// Row stride scalar multiply tests

	reset_values(size, a, b, c, c_copy);

	simple_smul(1, 1, &a[23], c, 4, 1);
	row_stride_smul_1x1(&a[23], c_copy, 4);
	test_buffer_diff(size, c, c_copy, "row_stride_smul_1x1");

	reset_values(size, a, b, c, c_copy);

	simple_smul(2, 1, &a[2], c, 2, 1);
	row_stride_smul_2x1(&a[2], c_copy, 2);
	test_buffer_diff(size, c, c_copy, "row_stride_smul_2x1");

	simple_smul(1, 2, &a[3], c, 6, 1);
	row_stride_smul_1x2(&a[3], c_copy, 6);
	test_buffer_diff(size, c, c_copy, "row_stride_smul_1x2");

	simple_smul(2, 2, &a[4], c, 6, 1);
	row_stride_smul_2x2(&a[4], c_copy, 6);
	test_buffer_diff(size, c, c_copy, "row_stride_smul_2x2");

	printf("\n\n---------- END ROW STRIDE TESTS ----------\n\n");
}