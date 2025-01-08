#!/bin/tcsh
#SBATCH --job-name=nextflow_test                              # job name
#SBATCH --partition=super                                     # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                             # number of nodes requested by user
#SBATCH --time=60-00:00:00                                    # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                             # standard output file name
#SBATCH --error=serialJob.%j.time                             # standard error output file name
#SBATCH --mail-user=Yihan.Wang@utsouthwestern.edu             # specify an email address
#SBATCH --mail-type=end                                       # send email when job status change (start, end, abortion and etc.)

module load nextflow/24.04.4
module load singularity/3.9.9

nextflow run -with-singularity pyspade_v0150.sif main.nf
