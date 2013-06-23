#!/bin/sh

#Gemm A/B/C N/T N/T
./dxter.x 5 C N N | tee BLASResults/resGemmNN
./dxter.x 5 C N T | tee BLASResults/resGemmNT
./dxter.x 5 C T N | tee BLASResults/resGemmTN
./dxter.x 5 C T T | tee BLASResults/resGemmTT

#Hemm
./dxter.x 17 C R L L | tee BLASResults/resHemmLL
./dxter.x 17 C R L U | tee BLASResults/resHemmLU
./dxter.x 17 C R R L | tee BLASResults/resHemmRL
./dxter.x 17 C R R U | tee BLASResults/resHemmRU

#Trsm 0/1/2 L/R L/U N/T
./dxter.x 3 2 L L N | tee BLASResults/resTrsmLLN
./dxter.x 3 2 L L T | tee BLASResults/resTrsmLLT
./dxter.x 3 2 L U N | tee BLASResults/resTrsmLUN
./dxter.x 3 2 L U T | tee BLASResults/resTrsmLUT
./dxter.x 3 2 R L N | tee BLASResults/resTrsmRLN
./dxter.x 3 2 R L T | tee BLASResults/resTrsmRLT
./dxter.x 3 2 R U N | tee BLASResults/resTrsmRUN
./dxter.x 3 2 R U T | tee BLASResults/resTrsmRUT

#Trmm A/C L/R L/U N/T
./dxter.x 18 C L L N | tee BLASResults/resTrmmLLN
./dxter.x 18 C L L T | tee BLASResults/resTrmmLLT
./dxter.x 18 C L U N | tee BLASResults/resTrmmLUN
./dxter.x 18 C L U T | tee BLASResults/resTrmmLUT
./dxter.x 18 C R L N | tee BLASResults/resTrmmRLN
./dxter.x 18 C R L T | tee BLASResults/resTrmmRLT
./dxter.x 18 C R U N | tee BLASResults/resTrmmRUN
./dxter.x 18 C R U T | tee BLASResults/resTrmmRUT

#Herk C/R L/U N/T/C
./dxter.x 15 R L N | tee BLASResults/resHerkLN
./dxter.x 15 R L T | tee BLASResults/resHerkLT
./dxter.x 15 R U N | tee BLASResults/resHerkUN
./dxter.x 15 R U T  | tee BLASResults/resHerkUT

#Her2k C/R L/U N/T/C
./dxter.x 16 R L N | tee BLASResults/resHer2kLN
./dxter.x 16 R L T | tee BLASResults/resHer2kLT
./dxter.x 16 R U N | tee BLASResults/resHer2kUN
./dxter.x 16 R U T  | tee BLASResults/resHer2kUT
