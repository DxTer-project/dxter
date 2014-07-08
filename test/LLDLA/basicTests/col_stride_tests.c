#include "col_stride_lldla_primitives.h"
#include "utils.h"

void column_stride_tests()	{
	int size = 50;
	double a[size];
	double b[size];
	double c[size];
	double c_copy[size];

	reset_values(size, a, b, c, c_copy);

	printf("\n\n---------- START COLUMN STRIDE TESTS ----------\n\n");

	test_buffer_diff(size, c, c_copy, "sanity check");

	// General stride add tests

	simple_add(2, 1, a, 1, 4, c, 1, 4);
       	col_stride_add_2x1(a, 4, c_copy, 4);
       	test_buffer_diff(size, c, c_copy, "col_stride_add_2x1 same row and col stride");

	reset_values(size, a, b, c, c_copy);

	simple_add(1, 2, a, 1, 2, c, 1, 6);
       	col_stride_add_1x2(a, 2, c_copy, 6);
       	test_buffer_diff(size, c, c_copy, "col_stride_add_1x2");

	reset_values(size, a, b, c, c_copy);

	simple_add(2, 2, a, 1, 2, c, 1, 6);
       	col_stride_add_2x2(a, 2, c_copy, 6);
       	test_buffer_diff(size, c, c_copy, "col_stride_add_2x2");

	reset_values(size, a, b, c, c_copy);

	double alpha = 3.2;
	simple_smul(2, 1, &alpha, c, 1, 6);
       	col_stride_smul_2x1(&alpha, c_copy, 6);
       	test_buffer_diff(size, c, c_copy, "col_stride_smul_2x1");

	reset_values(size, a, b, c, c_copy);

	alpha = -2.3;
	simple_smul(1, 2, &alpha, c, 1, 6);
       	col_stride_smul_1x2(&alpha, c_copy, 6);
       	test_buffer_diff(size, c, c_copy, "col_stride_smul_1x2");

	reset_values(size, a, b, c, c_copy);

	alpha = 98.18;
	simple_smul(2, 2, &alpha, c, 1, 8);
       	col_stride_smul_2x2(&alpha, c_copy, 8);
       	test_buffer_diff(size, c, c_copy, "col_stride_smul_2x2");

	printf("\n\n---------- END COLUMN STRIDE TESTS ----------\n\n");
}
