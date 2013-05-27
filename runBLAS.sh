#!/bin/sh

#Gemm A/B/C N/T N/T
./driver.x 5 C N N | tee BLASResults/resGemmNN
./driver.x 5 C N T | tee BLASResults/resGemmNT
./driver.x 5 C T N | tee BLASResults/resGemmTN
./driver.x 5 C T T | tee BLASResults/resGemmTT

#Hemm
./driver.x 17 C R L L | tee BLASResults/resHemmLL
./driver.x 17 C R L U | tee BLASResults/resHemmLU
./driver.x 17 C R R L | tee BLASResults/resHemmRL
./driver.x 17 C R R U | tee BLASResults/resHemmRU

#Trsm 0/1/2 L/R L/U N/T
./driver.x 3 2 L L N | tee BLASResults/resTrsmLLN
./driver.x 3 2 L L T | tee BLASResults/resTrsmLLT
./driver.x 3 2 L U N | tee BLASResults/resTrsmLUN
./driver.x 3 2 L U T | tee BLASResults/resTrsmLUT
./driver.x 3 2 R L N | tee BLASResults/resTrsmRLN
./driver.x 3 2 R L T | tee BLASResults/resTrsmRLT
./driver.x 3 2 R U N | tee BLASResults/resTrsmRUN
./driver.x 3 2 R U T | tee BLASResults/resTrsmRUT

#Trmm A/C L/R L/U N/T
./driver.x 18 C L L N | tee BLASResults/resTrmmLLN
./driver.x 18 C L L T | tee BLASResults/resTrmmLLT
./driver.x 18 C L U N | tee BLASResults/resTrmmLUN
./driver.x 18 C L U T | tee BLASResults/resTrmmLUT
./driver.x 18 C R L N | tee BLASResults/resTrmmRLN
./driver.x 18 C R L T | tee BLASResults/resTrmmRLT
./driver.x 18 C R U N | tee BLASResults/resTrmmRUN
./driver.x 18 C R U T | tee BLASResults/resTrmmRUT

#Herk C/R L/U N/T/C
./driver.x 15 R L N | tee BLASResults/resHerkLN
./driver.x 15 R L T | tee BLASResults/resHerkLT
./driver.x 15 R U N | tee BLASResults/resHerkUN
./driver.x 15 R U T  | tee BLASResults/resHerkUT

#Her2k C/R L/U N/T/C
./driver.x 16 R L N | tee BLASResults/resHer2kLN
./driver.x 16 R L T | tee BLASResults/resHer2kLT
./driver.x 16 R U N | tee BLASResults/resHer2kUN
./driver.x 16 R U T  | tee BLASResults/resHer2kUT
