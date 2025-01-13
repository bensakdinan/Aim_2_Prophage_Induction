#!/bin/sh
#SBATCH --mail-user=ben.sakdinan@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=propagate_with_VIBRANT
#SBATCH --output=%x-%j.out
#SBATCH --time=5:00:00
#SBATCH --tasks=6
#SBATCH --mem=100G
#SBATCH -o propagate_with_VIBRANT.o%j
#SBATCH --array=1-344

# sbatch --account=ctb-shapiro propagate_with_VIBRANT.sh

# Modules
module load StdEnv/2020
module load bowtie2 samtools
module load python

# SLURM array assignment
vibrant_output_dir=/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_vibrant/vibrant_output
prophage=$(cat /home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/samples.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID') 
echo $prophage

# Assigning Propagate input fields. 
scaffold=/home/bens/projects/ctb-shapiro/bens/prophage_induction/01_rerun_megahit_output_alinedata_01.10/${prophage}/final.contigs.fa
prophage_coordinates=$vibrant_output_dir/$prophage/VIBRANT_1kb_final.contigs/VIBRANT_results_1kb_final.contigs/VIBRANT_integrated_prophage_coordinates_1kb_final.contigs.tsv
reads=/home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/01_reads/${prophage}-QUALITY_PASSED

# Run Propagate

# Strict threshold
Propagate -f $scaffold -v $prophage_coordinates -r ${reads}_R1.fastq ${reads}_R2.fastq -t 12 -o /home/bens/projects/ctb-shapiro/bens/prophage_induction/04_propagate_with_VIBRANT/$prophage --clean

# Loose threshold
#Propagate -f $scaffold -v $prophage_coordinates -r ${reads}_R1.fastq ${reads}_R2.fastq -e 0.60 -c 1.50 -t 12 -o /home/bens/projects/ctb-shapiro/bens/prophage_induction/04_propagate_0.60+1.50_threshold/$prophage --clean

echo "Propagate done"
