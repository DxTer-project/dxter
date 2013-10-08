#!/bin/sh
export GOMP_CPU_AFFINITY="0 12 4 16 8 20 1 13 5 17 9 21 2 14 6 18 10 22 3 15 7 19 11 23"

# ./test_gemm.x N N | tee results/24/resGemmNN
# ./test_gemm.x T T | tee results/24/resGemmTT

# ./test_trmm.x L L N | tee results/24/resTrmmLLN
# ./test_trmm.x R L N | tee results/24/resTrmmRLN

# ./test_trsm.x L L N | tee results/24/resTrsmLLN
# ./test_trsm.x R L N | tee results/24/resTrsmRLN

# ./test_hemm.x L L | tee results/24/resHemmLL
# ./test_hemm.x R L | tee results/24/resHemmRL

./test_herk.x L N | tee results/24/resHerkLN
./test_herk.x U N | tee results/24/resHerkUN

./test_her2k.x L N | tee results/24/resHer2kLN
./test_her2k.x U N | tee results/24/resHer2kUN
