#!/bin/sh
#SBATCH --mail-user=ben.sakdinan@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=genomad_with_array
#SBATCH --output=%x-%j.out
#SBATCH --time=5:00:00
#SBATCH --tasks=20
#SBATCH --mem=100G
#SBATCH -o genomad_with_array.o%j
#SBATCH --array=1-345

# sbatch --account=ctb-shapiro geNomad_with_array.sh

working_dir=/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad

# Modules and environments
module load python/3.10 mmseqs2 aragorn

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip

# SLURM array 
metagenomes_directory=/home/bens/projects/ctb-shapiro/bens/prophage_induction/01_rerun_megahit_output_alinedata_01.10
metagenome=$(cat /home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/samples.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID') 
echo $metagenome

# Make a subdirectory for the current Vc genome
# if mkdir error, add '-p'
if [ -d "$working_dir/$metagenome" ]; then
    rm -rf "$working_dir/$metagenome"
fi
mkdir -p "$working_dir/$metagenome"

# Genomad input fields
database=/home/bens/projects/ctb-shapiro/bens/software/genomad/genomad_db
out=$working_dir/$metagenome

# Run genomad  
genomad end-to-end  --cleanup --splits 8 --enable-score-calibration $metagenomes_directory/$metagenome/final.contigs.fa $out $database

deactivate
