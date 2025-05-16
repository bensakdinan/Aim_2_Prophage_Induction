# Prophage induction re-analysis

## CheckV
Run checkV to identify completeness of each prophage contig. For each VC that was previously constructed, we only want to keep one representative contig per cluster. The contig with the highest completeness score will be kept.

```
#!/usr/bin/bash
#SBATCH -N 1                           
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=8:30:00 
#SBATCH -o checkv.o%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ben.sakdinan@mail.mcgill.ca

# sbatch --account=ctb-shapiro checkv.sh

module load StdEnv/2020 python/3.11.5 diamond/2.0.4 hmmer/3.3.2 prodigal/2.6.3 prodigal-gv/2.6.3 

source /home/bens/projects/ctb-shapiro/bens/software/checkv/checkv_env/bin/activate

export CHECKVDB=/project/def-corinnem/databases/checkv-db-v1.4

input_file=/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad/all_prophages_seqs.faa
output_dir=/home/bens/projects/ctb-shapiro/bens/prophage_induction/08_checkv
checkv end_to_end $input_file $output_dir -t 16

deactivate 

```
