# vConTACT2 analysis in R

## Import libraries
```r
library(dplyr)
library(ggplot2)
```

## Filter the dataset to only keep my contigs, not the reference DB contigs
```r
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

## Distribution of cluster sizes
```r
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

## Per cluster, how many distinct samples are present 
```r
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
