#!/bin/sh
export GOMP_CPU_AFFINITY="0 12 4 16 8 20 1 13 5 17 9 21 2 14 6 18 10 22 3 15 7 19 11 23"

./test_gemm.x N N | tee results/16/resGemmNN
./test_trmm.x L L N | tee results/16/resTrmmLLN
./test_herk.x L N | tee results/16/resHerkLN