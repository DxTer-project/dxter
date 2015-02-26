#!/bin/bash
#SBATCH -J twoSidedTrmm           # job name
#SBATCH -o twoSidedTrmm.o%j       # output and error file name (%j expands to jobID)
#SBATCH -e twoSidedTrmm.e%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1              # total number of mpi tasks requested
#SBATCH -p normal     # queue (partition) -- normal, development, etc.
#SBATCH -t 02:30:00        # run time (hh:mm:ss)
#SBATCH --mail-user=bamarker@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

#./test_twoSidedTrmm.x | tee outputTwoSidedTrmm
./test_twoSidedTrmm_NoBLIS.x | tee outputTwoSidedTrmmNoBLIS