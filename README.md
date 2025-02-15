# Aim_1_Prophage_Induction
This markdown contains my step by step bioinformatic protocol

## Dataset

The dataset I am working with is a previously established dataset from our lab. Illumina short reads were acquired and stored at 
`/home/bens/projects/ctb-shapiro/bens/aline_data_07_02_24/01_data/01_reads`


## Metagenome Assembly

Short reads were assembled using MEGAHIT. 
https://github.com/voutcn/megahit

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
https://github.com/apcamargo/genomad

GeNomad outputs can be found at `/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad`
Prophage data can be found at `./${Sample}/final.contigs_find_proviruses`
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

#### Prophage Coordinates
A later tool in my bioinformatic pipeline, PropagAtE, requires a prophage coordinates file. This is a .tsv file that contains the scaffold name that contains the prophage, the name of the specific propahge fragment name, the start nucleotide, and the stop nucleotide. GeNomad has this information within a larger file named final.contigs_provirus.tsv, but PropagAtE needs a .tsv file with only these data. To make this file, the code below was used:
```bash
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

```

### VIBRANT
VIBRANT is another viral contig identification tool. 
https://github.com/AnantharamanLab/VIBRANT 

This tool was run by a colleague in the lab, Dr. Anshul Sinha, the output can be found at `/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_vibrant/vibrant_output`

## Detection of Prophage Induction - PropagAtE
PropagAtE is a tool to identify actively replicating prophages. Propagate compares the depth of prophage reads to entire host region reads to classify a phage as active or inactive at a certain read depth ratio. These thresholds can be manually changed by the user.
https://github.com/AnantharamanLab/PropagAtE

Depending on whether GeNomad or VIBRANT was used, a slightly different PropagAtE script must be used. I ran PropagAtE at different a conservative (default) threshold and liberal threshold by adjusting these flags `-e 0.60 -c 1.50`

GeNomad conservative threshold output: `/home/bens/projects/ctb-shapiro/bens/prophage_induction/04_propagate`
GeNomad liberal threshold output: `/home/bens/projects/ctb-shapiro/bens/prophage_induction/04_propagate_0.60+1.50_threshold`
VIBRANT default threshold output: `/home/bens/projects/ctb-shapiro/bens/prophage_induction/04_propagate_with_VIBRANT`

For GeNomad:
```bash
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
```

For VIBRANT:
```bash
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
```

#### Prophage:host read ratio
PropagAtE will classify which prophages are actively replicating, however, I want to also analyze the raw data and see the exact ratios of prophage:host read ratios, not just the classifications of "active" or "inactive". This script below takes a detailed PropagAtE results output and organizes a .csv file containing each sample's prophage:host read ratio for each prophage in that given sample.

Prophage:host read ratios can be found at: `/home/bens/projects/ctb-shapiro/bens/prophage_induction/05_analysis`
```bash
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
```

## VContact2: Clustering by proteins

#### Concatenating prophage proteins from geNomad output
```bash
# Get .faa protein files from geNomad outputs and concatenate into a single .faa file
# geNomad outputs at /home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad/${SAMPLE}
for sample_dir in "/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad"/*/; do
    protein_file="${sample_dir}final.contigs_find_proviruses/final.contigs_provirus_proteins.faa"

    if [[ -f "$protein_file" ]]; then
        cat "$protein_file" >> "all_unique_proteins.faa"
    else
        echo "Warning: File not found - $protein_file"
    fi
done

# Remove trailing strings from protein identifier lines
awk '{if ($0 ~ /^>/) {print $1} else {print $0}}' all_unique_proteins.faa > concatenated_prophage_proteins.faa
```

#### Making gene-to-genome mapping file
```bash
output_csv="test_g2g_vcontact_00.csv"

# Write the headers row to output file
echo "protein_id,contig_id,keywords" > "$output_csv"

# Loop through each sample directory
for sample_dir in ./*/; do
    # Extract the sample name from the directory path
    sample_name=$(basename "$sample_dir")

    # Path to the geNomad provirus genes .tsv file
    genes_file="${sample_dir}final.contigs_find_proviruses/final.contigs_provirus_genes.tsv"

    # Check if the file exists
    if [[ -f "$genes_file" ]]; then
        echo "Processing $genes_file..."

        # Extract the required columns and modify contig_id
        awk -F'\t' -v sample="$sample_name" 'NR > 1 {
            # Full protein ID with sample name
            protein_id = sample "_" $1;

            # Contig ID is the protein ID without sample name and trailing _# 
            contig_id = $1;  
            sub(/_[0-9]+$/, "", contig_id)  # Remove the final underscore + number

            # Annotation description (protein function info)
            split($20, keywords, ",")
            first_keyword = keywords[1];

            print protein_id "," contig_id "," first_keyword
        }' "$genes_file" >> "$output_csv"
    else
        echo "Warning: File not found in $genes_file"
    fi
done

```

#### Running vConTACT2
First create virtualenv for vConTACT2
- scipy-stack/2020a is necessary to load pandas/0.24.2, otherwise pandas dependency error will occur
```bash
module load StdEnv/2020 gcc python/3.8.10 scipy-stack/2020a hdf5/1.12.1 diamond blast+ mcl java
virtualenv --no-download --clear ~/ENV && source ~/ENV/bin/activate
pip install 'numpy<1.24' git+https://bitbucket.org/MAVERICLab/vcontact2.git
vcontact2 -h
wget https://github.com/ntamas/cl1/releases/download/1.2/cluster_one-1.2.jar -O $VIRTUAL_ENV/bin/cluster_one-1.2.jar
```
Run vConTACT2 script
```bash
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

module load StdEnv/2020 gcc python/3.8.10 scipy-stack/2020a hdf5/1.12.1 diamond blast+ mcl java 

source ~/ENV/bin/activate

# define paths
input_protein_file=/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad/concatenated_prophage_proteins.faa
input_gene_to_genome=/home/bens/projects/ctb-shapiro/bens/prophage_induction/03_genomad/g2g_vcontact.csv
output_dir=/home/bens/projects/ctb-shapiro/bens/prophage_induction/06_vcontact2

# Run vcontact2
vcontact2 --raw-proteins $input_protein_file --proteins-fp $input_gene_to_genome --db 'ProkaryoticViralRefSeq94-Merged' --output-dir $output_dir --c1-bin /home/bens/ENV/bin/cluster_one-1.2.jar

deactivate
```

### Analyze vConTACT2 output in R

- First map each contig to its original sample (in bash)
```bash
# Map each contig to its original sample

echo "contig,Sample" > "sample_by_contig.csv"

for sample_dir in */; docd ./
    sample_name="${sample_dir%/}"  # Remove trailing slash to get sample name
    fasta_file="${sample_dir}final.contigs_find_proviruses/final.contigs_provirus.fna"

    if [[ -f "$fasta_file" ]]; then
        # Extract contig names and write to CSV
        grep "^>" "$fasta_file" | sed 's/^>//' | awk -v sample="$sample_name" '{print $1 "," sample}' >> "sample_by_contig.csv"
    fi
done
```

#### Import libraries
```{r}
library(dplyr)
library(ggplot2)
```

#### Filter the dataset to only keep my contigs, not the reference DB contigs
```{r}
filtered_genome_by_genome <- genome_by_genome_overview %>% 
  filter(grepl("^k", Genome))

# Total number of contigs 
nrow(filtered_genome_by_genome)

# Filter out singletons and outliers 
clustered_genome_by_genome <- filtered_genome_by_genome %>%
  filter(`VC Status` == "Clustered")

# Total number and proportion of clustered contigs
nrow(clustered_genome_by_genome)
nrow(clustered_genome_by_genome) / nrow(filtered_genome_by_genome)
```

- 2941/4055 (72.53%) of contigs are clustered (not singletons or outliers)

#### Distribution of cluster sizes
```{r}
# How many unique clusters exist
clusters <- unique(clustered_genome_by_genome$VC)
length(clusters)

# Get each unique VC and its VC size
vclusters <- clustered_genome_by_genome %>%
  select(VC, `VC Size`) %>%
  distinct()

# Plot histogram by largest > smallest
vclusters <- vclusters %>%
  arrange(desc(`VC Size`))

vclusters$VC <- factor(vclusters$VC, levels = vclusters$VC)

# Plot
ggplot(vclusters, aes(x = VC, y = `VC Size`)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Distribution of VCs by VC Size", x = "VC (477 Total)", y = "VC Size") +
  theme(axis.text.x = element_blank())

# How many are doublets or triplets
print("Number of doublets")
sum(vclusters$`VC Size` == 2)

print("Number of triplets")
sum(vclusters$`VC Size` == 3)
```

#### Per cluster, how many distinct samples are present 
```{r}
# Adding sample ID and Dehydration_Status per contig
samples_in_vclusters <- sample_by_contig %>%
  left_join(X00_metadata %>%
              select(Sample, Dehydration_Status), by = "Sample")

# Adding cluster ID  ("Genome") per contig
samples_in_vclusters <- samples_in_vclusters %>%
  left_join(filtered_genome_by_genome %>% select(Genome, VC), by = "Genome")

# Filter out rows with NA values
samples_in_vclusters <- samples_in_vclusters %>% # Removes unclustered contigs
  filter(!is.na(VC) )
samples_in_vclusters <- samples_in_vclusters %>% # Removes positive control 
  filter(!is.na(Dehydration_Status) )

# Remaining total number of contigs after clustering and filtering
nrow(samples_in_vclusters)

# For each VC cluster, how many contigs are present at each Dehydration Status
VC_vs_numSamples <- samples_in_vclusters %>%
  group_by(VC) %>%
  summarise(
    Mild = sum(Dehydration_Status == 1, na.rm = TRUE),
    Moderate = sum(Dehydration_Status == 2, na.rm = TRUE),
    Severe = sum(Dehydration_Status == 3, na.rm = TRUE)
  ) %>% 
  ungroup()
```

