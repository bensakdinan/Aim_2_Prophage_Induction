#!/bin/sh
#SBATCH --mail-user=ben.sakdinan@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=prophage_coordinates
#SBATCH --output=%x-%j.out
#SBATCH --time=1:00
#SBATCH -o prophage_coordinates.o%j
#SBATCH --array=1-345

# sbatch --account=ctb-shapiro prophage_coordinates.sh


prophage=$(cat /home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/samples.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID') 
genomad_output_dir=/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad
echo $prophage

prophage_dir=$genomad_output_dir/$prophage/final.contigs_find_proviruses
file=$prophage_dir/final.contigs_provirus.tsv

cd $prophage_dir
if [ -f "${prophage}_prophage_coordinates.tsv" ]; then
    rm ${prophage}_prophage_coordinates.tsv
fi 

# Prophage coordinates for MEGAHIT assemblies
# I need to write code that will need to go to 01_megahit folder and get the name of the entire header line AND THEN assign it to a variable, VARIABLE_OF_HEADER, which will be awk printed in the below block of code
if [ -f "$file" ]; then 
    echo -e "scaffold\tfragment\tstart\tstop" > ${prophage}_prophage_coordinates.tsv

    metagenome=/home/bens/projects/ctb-shapiro/bens/prophage_induction/01_rerun_megahit_output_alinedata_01.10/${prophage}/final.contigs.fa

    awk 'NR > 1 {print $2}' "$file" | while read -r contig; do
        scaffold_header=$(grep -m 1 "$contig" "$metagenome" | sed 's/^>//')

        awk -v header="$scaffold_header" -v contig="$contig" 'NR > 1 && $2 == contig {
            print header "\t" $1 "\t" $3 "\t" $4
        }' "$file" >> ${prophage}_prophage_coordinates.tsv
    done
fi
