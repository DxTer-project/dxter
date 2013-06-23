#!/bin/bash
#SBATCH -J twoSidedTrsm           # job name
#SBATCH -o twoSidedTrsm.o%j       # output and error file name (%j expands to jobID)
#SBATCH -e twoSidedTrsm.e%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1              # total number of mpi tasks requested
#SBATCH -p normal     # queue (partition) -- normal, development, etc.
#SBATCH -t 02:30:00        # run time (hh:mm:ss)
#SBATCH --mail-user=bamarker@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

#./test_twoSidedTrsm.x | tee outputTwoSidedTrsm

./test_twoSidedTrsm_NoBLIS.x | tee outputTwoSidedTrsmNoBLIS