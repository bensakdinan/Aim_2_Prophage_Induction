# Prophage induction re-analysis

## CheckV
Run checkV to identify completeness of each prophage contig. For each VC that was previously constructed, we only want to keep one representative contig per cluster. The contig with the highest completeness score will be kept.

```bash
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

## Concatenate prophage coordiantes for these representative contigs 
```bash
#!/bin/bash

cd /home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad
BASE_DIR="/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad"
OUTPUT_FILE="representative_prophage_coordinates.tsv"

# Initialize output file with header (copied from first file)
FIRST_SAMPLE=$(awk -F, 'NR==2 {print $2}' /home/bens/projects/ctb-shapiro/bens/prophage_induction/08_checkv/representative_prophage_contigs.csv)
FIRST_COORD_FILE="$BASE_DIR/${FIRST_SAMPLE}/final.contigs_find_proviruses/${FIRST_SAMPLE}_prophage_coordinates.tsv"
echo -e "scaffold\tfragment\tstart\tstop" > "$OUTPUT_FILE"
head -n 1 "$FIRST_COORD_FILE" > "$OUTPUT_FILE"

# Iterate through the CSV file, skipping header
tail -n +2 /home/bens/projects/ctb-shapiro/bens/prophage_induction/08_checkv/representative_prophage_contigs.csv | while IFS=, read -r _ Genome _ Sample _; do
    # Remove surrounding quotation marks from values
    Genome=$(echo "$Genome" | tr -d '"')
    Sample=$(echo "$Sample" | tr -d '"')

    COORD_FILE="$BASE_DIR/${Sample}/final.contigs_find_proviruses/${Sample}_prophage_coordinates.tsv"
    echo $COORD_FILE

    if [[ -f "$COORD_FILE" ]]; then
        awk -v contig="$Genome" -F'\t' '$2 == contig' "$COORD_FILE" >> "$OUTPUT_FILE"
    else
        echo "Warning: File not found: $COORD_FILE" >&2
    fi
done

```
