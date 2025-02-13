#!/bin/sh
#SBATCH --mail-user=ben.sakdinan@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=vcontact2
#SBATCH --output=%x-%j.out
#SBATCH --time=5:00:00
#SBATCH --tasks=6
#SBATCH --mem=64G
#SBATCH -o vcontact2.o%j

# sbatch --account=ctb-shapiro vcontact2.sh

# module load StdEnv/2020 gcc python/3.10 hdf5/1.12.1 diamond blast+ mcl java
module load StdEnv/2020 gcc python/3.8.10 scipy-stack/2020a hdf5/1.12.1 diamond blast+ mcl java 

source ~/ENV/bin/activate

# define paths
input_protein_file=/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad/concatenated_prophage_proteins.faa
input_gene_to_genome=/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad/g2g_vcontact.csv
output_dir=/home/bens/projects/ctb-shapiro/bens/prophage_induction/06_vcontact2_rerun1

# Run vcontact2
vcontact2 --raw-proteins $input_protein_file --proteins-fp $input_gene_to_genome --db 'ProkaryoticViralRefSeq94-Merged' --output-dir $output_dir --c1-bin /home/bens/ENV/bin/cluster_one-1.2.jar

deactivate
