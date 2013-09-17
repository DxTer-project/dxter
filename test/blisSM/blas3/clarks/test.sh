#!/bin/sh
export GOMP_CPU_AFFINITY="0 12 4 16 8 20 1 13 5 17 9 21 2 14 6 18 10 22 3 15 7 19 11 23"

./test_gemm.x N N | tee results/16/resGemmNN
./test_gemm.x T T | tee results/16/resGemmTT

./test_trmm.x L L N | tee results/16/resTrmmLLN
./test_trmm.x R L N | tee results/16/resTrmmRLN

./test_trsm.x L L N | tee results/16/resTrsmLLN
./test_trsm.x R L N | tee results/16/resTrsmRLN


./test_hemm.x L L | tee results/16/resHemmLL
./test_hemm.x R L | tee results/16/resHemmRL
