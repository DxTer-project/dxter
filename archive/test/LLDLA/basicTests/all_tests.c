#include "col_stride_tests.h"
#include "gen_stride_tests.h"
#include "row_stride_tests.h"
#include "utils.h"

int main()	{
  column_stride_tests();
  general_stride_tests();
  row_stride_tests();
  return 0;
}
