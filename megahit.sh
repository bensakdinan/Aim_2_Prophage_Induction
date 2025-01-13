#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --tasks=12
#SBATCH --mem=100G
#SBATCH -o megahit.o%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ben.sakdinan@mail.mcgill.ca
#SBATCH --array=1-345

#How to run
# sbatch --account=def-shapiro-ab megahit.sh 

module load StdEnv/2023
module load megahit/1.2.9

PATIENT=$(cat /home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/samples.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID')

#Make output directory
#cd /home/bens/projects/ctb-shapiro/bens/prophage_induction/megahit_output_alinedata_01.10
#mkdir -p $PATIENT

#Run megahit 
READS=/home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/01_reads

megahit -1 $READS/${PATIENT}-QUALITY_PASSED_R1.fastq -2 $READS/${PATIENT}-QUALITY_PASSED_R2.fastq -t 16 -o /home/bens/projects/ctb-shapiro/bens/prophage_induction/rerun_megahit_output_alinedata_01.10/$PATIENT
