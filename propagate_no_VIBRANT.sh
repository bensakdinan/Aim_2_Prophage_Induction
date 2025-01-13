#!/bin/sh
#SBATCH --mail-user=ben.sakdinan@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=propagate_no_VIBRANT
#SBATCH --output=%x-%j.out
#SBATCH --time=5:00:00
#SBATCH --tasks=6
#SBATCH --mem=100G
#SBATCH -o propagate_no_VIBRANT.o%j
#SBATCH --array=1-6

# sbatch --account=ctb-shapiro propagate_no_VIBRANT.sh

# Modules
module load StdEnv/2020
module load bowtie2 samtools
module load python

# SLURM array assignment
genomad_output_dir=/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad
prophage=$(cat /home/bens/projects/ctb-shapiro/bens/scripts/misc_files/propagate_reruns.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID') 
echo $prophage

# Prophage output directory from geNomad
prophage_dir=$genomad_output_dir/$prophage/final.contigs_find_proviruses

# Assigning Propagate input fields. 
scaffold=/home/bens/projects/ctb-shapiro/bens/prophage_induction/01_rerun_megahit_output_alinedata_01.10/${prophage}/final.contigs.fa
prophage_coordinates=$prophage_dir/${prophage}_prophage_coordinates.tsv
reads=/home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/01_reads/${prophage}-QUALITY_PASSED

# Run Propagate
Propagate -f $scaffold -v $prophage_coordinates -r ${reads}_R1.fastq ${reads}_R2.fastq -e 0.60 -c 1.50 -t 12 -o /home/bens/projects/ctb-shapiro/bens/prophage_induction/04_propagate_0.60+1.50_threshold/$prophage --clean

echo "Propagate done"

#Propagate -f /home/bens/projects/ctb-shapiro/bens/EN23_VC_genomes_fasta/1376450.fasta -v /home/bens/scratch/geNomad_output_21.10/1376450/1376450_find_proviruses/1376450_prophage_coordinates.tsv -r /home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/01_reads/1376450-QUALITY_PASSED_R1.fastq /home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/01_reads/1376450-QUALITY_PASSED_R2.fastq -t 12 -o /home/bens/scratch/chuhan_propagate/tempdir_bugfixing --clean
