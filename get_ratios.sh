#!/bin/sh
#SBATCH --mail-user=ben.sakdinan@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=get_ratios
#SBATCH --output=%x-%j.out
#SBATCH --time=00:15:00
#SBATCH --tasks=6
#SBATCH --mem=16G
#SBATCH -o get_ratios.o%j
#SBATCH --array=1-344

# sbatch --account=ctb-shapiro get_ratios.sh

metagenome=$(cat /home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/samples.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID') 
echo $metagenome

log_file=/home/bens/projects/ctb-shapiro/bens/prophage_induction/04_propagate/${metagenome}/${metagenome}.tsv

cd /home/bens/projects/ctb-shapiro/bens/prophage_induction/05_analysis

awk -v metagenome="${metagenome}" 'BEGIN { OFS="\t" } NR > 1 { print metagenome, $8 }' "/home/bens/projects/ctb-shapiro/bens/prophage_induction/04_propagate/${metagenome}/${metagenome}.tsv" >> prophage-host_ratios.tsv
