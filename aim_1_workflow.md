# This markdown contains my step by step bioinformatic protocol

## Dataset

The dataset I am working with is a previously established dataset from our lab. Illumina short reads were acquired and stored at 
`/home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/01_reads`


## Metagenome Assembly

Short reads were assembled using MEGAHIT. 
```bash
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
```
Assembled metagenomes for each sample were stored in their own subdirectories at `/home/bens/projects/ctb-shapiro/bens/prophage_induction/01_rerun_megahit_output_alinedata_01.10`. 
Contigs assembled in each genome were stored in the file `final.contigs.fa` in each sample's subdirectory 

### Contigs >=1kb filtering 
Only contigs >=1kb were kept for further analysis. For some further analyses, these contigs could be pooled into one .fasta file stored at `reran_concatenated_contigs_1000kb_11.10.fasta.gz`
```bash
module load seqkit/2.5.1

seqkit seq -m 1000 -g /home/bens/projects/ctb-shapiro/bens/prophage_induction/rerun_megahit_output_alinedata_01.10/new_concatenated_contigs_11.10.fasta -o /home/bens/projects/ctb-shapiro/bens/prophage_induction/rerun_megahit_output_alinedata_01.10/reran_concatenated_contigs_1000kb_11.10.fasta
```
## Identification of viral contigs and prophages 

### GeNomad 
GeNomad is a tool to identify viral contigs. An functionality of this tool also identifies prophage contigs. 
```bash
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
```



