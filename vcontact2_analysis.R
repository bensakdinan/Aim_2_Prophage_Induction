---
title: "vcontact2_analysis"
author: "Ben Sakdinan"
date: "2025-02-14"
output: html_document
---

vConTACT2 analysis in R

## Import libraries
```{r}
library(dplyr)
library(ggplot2)
```


## Data setup

### Filter the dataset to only keep my contigs, not the reference DB contigs
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

### Distribution of cluster sizes
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
print("Number of doublets")
sum(vclusters$`VC Size` == 4)
print("Number of triplets")
sum(vclusters$`VC Size` == 5)
print("Number of triplets")
sum(vclusters$`VC Size` > 5)
```

### Per cluster, how many distinct samples are present 
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

VC_vs_numSamples <- merge(VC_vs_numSamples, vclusters, by = 'VC')

# unstratified_clusters <- sum(rowSums(VC_vs_numSamples[, c("Mild", "Moderate", "Severe")] > 0) == 1)
# print(unstratified_clusters)
```

### What is the distrubution of cluster sizes among these 292 stratified clusters
```{r}
multi_category_VCs <- VC_vs_numSamples[rowSums(VC_vs_numSamples[, c("Mild", "Moderate", "Severe")] > 0) > 1, ]
nrow(multi_category_VCs)

# Plot histogram of cluster size distribution
# Get each unique VC and its VC size
stratified_clusters <- multi_category_VCs %>%
  select(VC, `VC Size`) %>%
  distinct()

# Plot histogram by largest > smallest
stratified_clusters <- stratified_clusters %>%
  arrange(desc(`VC Size`))

stratified_clusters$VC <- factor(stratified_clusters$VC, levels =
                                   stratified_clusters$VC)

# Plot 
ggplot(stratified_clusters, aes(x = VC, y = `VC Size`)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Distribution of VCs by VC Size", x = "VC", y = "VC Size") +
  theme(axis.text.x = element_blank())

# Distrbution of cluster sizes
print("Number of doublets")
sum(stratified_clusters$`VC Size` == 2)
print("Number of triplets")
sum(stratified_clusters$`VC Size` == 3)
print("Number of quartets")
sum(stratified_clusters$`VC Size` == 4)
print("Number of quintets")
sum(stratified_clusters$`VC Size` == 5)
print("Number of clusters >6")
sum(stratified_clusters$`VC Size` > 5)
sum(stratified_clusters$`VC Size` > 10)
```



## Dehydration Status
#### Assess Prophage induction in doublets (and binned by Mild/Non-Mild or Non-Severe/Severe)

Doublet mild/nonmild
```{r}
# Make new dataframe of doublets across samples only and bin by mild/non-mild
VC_vs_numSamples_tmp <- VC_vs_numSamples %>% select(-`VC Size`)
classification_target <- merge(VC_vs_numSamples_tmp, stratified_clusters, by = 'VC')
classification_target <- classification_target %>%
  mutate(NonMild = Moderate + Severe)
classification_target <- classification_target %>%
  filter((Mild+NonMild) == 2 & Mild != 0  & NonMild != 0)
```

```{r}
# In shell, make propagate_by_contigs

# Merge data from propagate results into sample_by_contig
# propagate_by_contigs <- merge(sample_by_contig, propagate_results_by_prophage_conservative, by = 'Genome')
propagate_by_contigs <- merge(sample_by_contig, propagate_results_by_prophage_liberal, by = 'Genome')

# Merge sample dehydration status to filtered_genome_by_genome
sample_and_status <- X00_metadata %>%
  select(Sample, Dehydration_Status)

# Merge filtered genome with sample_and_status, then remove duplicate columns, active and prophage-host_ratio

filtered_genome_by_genome <- genome_by_genome_overview %>% 
  filter(grepl("^k", Genome))
filtered_genome_by_genome <- merge(filtered_genome_by_genome, propagate_by_contigs, by = 'Genome')
filtered_genome_by_genome <- merge(filtered_genome_by_genome, sample_and_status, by = 'Sample')
filtered_genome_by_genome <- filtered_genome_by_genome %>% select(-active, -`prophage-host_ratio`)

final_contig_data <- merge(filtered_genome_by_genome, propagate_results_by_prophage_liberal, by = 'Genome')
final_contig_data <- final_contig_data %>%
  select(
    Genome, 
    VC,
    Sample, 
    Dehydration_Status,
    active, 
    `prophage-host_ratio`
  )
final_contig_data <- final_contig_data %>%
  filter(!is.na(VC))

clusters_and_contigs <- data.frame(
  VC = classification_target$VC,
  mild_status = rep(NA_character_, nrow(classification_target)),
  mild_ratio = rep(NA_real_, nrow(classification_target)),
  nonmild_status = rep(NA_character_, nrow(classification_target)),
  nonmild_ratio = rep(NA_real_, nrow(classification_target)),
  stringsAsFactors = FALSE
  )

for (i in seq_len(nrow(clusters_and_contigs))) {
  current_cluster <- clusters_and_contigs$VC[i]
  
  for (j in seq_len(nrow(final_contig_data))) {
    # print(j)
    if (final_contig_data$VC[j] == current_cluster) {
    current_row <- final_contig_data[j, ]
    
      if (current_row$Dehydration_Status == 1) {
        clusters_and_contigs$mild_status[i] <- current_row$active
        clusters_and_contigs$mild_ratio[i]  <- current_row$`prophage-host_ratio`
      } else {
        clusters_and_contigs$nonmild_status[i] <- current_row$active
        clusters_and_contigs$nonmild_ratio[i]  <- current_row$`prophage-host_ratio`
      }
    }
  }
}

```

#### Non-Severe/Severe
```{r}
# Make new dataframe of doublets across samples only and bin by mild/non-mild
VC_vs_numSamples_tmp <- VC_vs_numSamples %>% select(-`VC Size`)
classification_target <- merge(VC_vs_numSamples_tmp, stratified_clusters, by = 'VC')
classification_target <- classification_target %>%
  mutate(NonSevere = Moderate + Mild)
classification_target <- classification_target %>%
  filter((Severe+NonSevere) == 2 & Severe != 0  & NonSevere != 0)
```

```{r}
# In shell, make propagate_by_contigs

# Merge data from propagate results into sample_by_contig
# propagate_by_contigs <- merge(sample_by_contig, propagate_results_by_prophage_conservative, by = 'Genome')
propagate_by_contigs <- merge(sample_by_contig, propagate_results_by_prophage_liberal, by = 'Genome')

# Merge sample dehydration status to filtered_genome_by_genome
sample_and_status <- X00_metadata %>%
  select(Sample, Dehydration_Status)

# Merge filtered genome with sample_and_status, then remove duplicate columns, active and prophage-host_ratio

filtered_genome_by_genome <- genome_by_genome_overview %>% 
  filter(grepl("^k", Genome))
filtered_genome_by_genome <- merge(filtered_genome_by_genome, propagate_by_contigs, by = 'Genome')
filtered_genome_by_genome <- merge(filtered_genome_by_genome, sample_and_status, by = 'Sample')
filtered_genome_by_genome <- filtered_genome_by_genome %>% select(-active, -`prophage-host_ratio`)

final_contig_data <- merge(filtered_genome_by_genome, propagate_results_by_prophage_liberal, by = 'Genome')
final_contig_data <- final_contig_data %>%
  select(
    Genome, 
    VC,
    Sample, 
    Dehydration_Status,
    active, 
    `prophage-host_ratio`
  )
final_contig_data <- final_contig_data %>%
  filter(!is.na(VC))

clusters_and_contigs <- data.frame(
  VC = classification_target$VC,
  severe_status = rep(NA_character_, nrow(classification_target)),
  severe_ratio = rep(NA_real_, nrow(classification_target)),
  nonsevere_status = rep(NA_character_, nrow(classification_target)),
  nonsevere_ratio = rep(NA_real_, nrow(classification_target)),
  stringsAsFactors = FALSE
  )

for (i in seq_len(nrow(clusters_and_contigs))) {
  current_cluster <- clusters_and_contigs$VC[i]
  
  for (j in seq_len(nrow(final_contig_data))) {
    # print(j)
    if (final_contig_data$VC[j] == current_cluster) {
    current_row <- final_contig_data[j, ]
    
      if (current_row$Dehydration_Status == 3) {
        clusters_and_contigs$severe_status[i] <- current_row$active
        clusters_and_contigs$severe_ratio[i]  <- current_row$`prophage-host_ratio`
      } else {
        clusters_and_contigs$nonsevere_status[i] <- current_row$active
        clusters_and_contigs$nonsevere_ratio[i]  <- current_row$`prophage-host_ratio`
      }
    }
  }
}

clusters_and_contigs <- clusters_and_contigs %>%
  filter(!is.na(severe_ratio))
clusters_and_contigs <- clusters_and_contigs %>%
  filter(!is.na(nonsevere_ratio))

```

```{r}
clusters_and_contigs$difference <- 
  (clusters_and_contigs$severe_ratio) - (clusters_and_contigs$nonsevere_ratio)
clusters_and_contigs <- clusters_and_contigs %>%
  arrange(desc(difference))

ggplot(clusters_and_contigs, aes(x = reorder(VC, -difference), y = difference, fill = difference > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue")) +
  labs(x = "VC", y = "Difference", title = "Difference in prophage-host read ratio, Severe - NonSevere (per VC)")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2))

# Test normality -> Not normally distributed
shapiro.test(clusters_and_contigs$difference)

# Stats testing
wilcox.test(clusters_and_contigs$difference, mu = 0, alternative = "two.sided")
```

#### Assess propahge induduction across all samples, averaged if theyre greater than a doublet

Mild vs NonMild
```{r}
# Make df to store all VCs (that are present in at least two levels of dehydration)
VC_vs_numSamples_tmp <- VC_vs_numSamples %>% select(-`VC Size`)
VC_vs_numSamples_tmp <- VC_vs_numSamples_tmp %>%
  mutate(NonMild = Moderate + Severe)
VC_vs_numSamples_tmp <- VC_vs_numSamples_tmp %>%
  mutate(NonSevere = Moderate + Mild)
classification_target <- merge(VC_vs_numSamples_tmp, stratified_clusters, by = 'VC')


# Initializing columns with NA values
# classification_target$mild_contigs <- vector("list", nrow(classification_target))
# classification_target$mild_ratios <- vector("list", nrow(classification_target))
# classification_target$mild_active <- vector("list", nrow(classification_target))

# classification_target$nonmild_contigs <- vector("list", nrow(classification_target))
# classification_target$nonmild_ratios <- vector("list", nrow(classification_target))
# classification_target$nonmild_active <- vector("list", nrow(classification_target))


# Then need to initialize columns for the average value
classification_target$avg_mild_ratios    <- rep(NA, nrow(classification_target))
classification_target$avg_mild_active    <- rep(NA, nrow(classification_target))
classification_target$avg_nonmild_ratios <- rep(NA, nrow(classification_target))
classification_target$avg_nonmild_active <- rep(NA, nrow(classification_target))

for (i in seq_len(nrow(classification_target))) {
  current_cluster <- as.character(classification_target$VC[i])
  
  mild_contigs <- character(0)
  mild_ratios <- numeric(0)
  mild_active <- character(0)
  
  nonmild_contigs <- character(0)
  nonmild_ratios <- numeric(0)
  nonmild_active <- character(0)

  for (j in seq_len(nrow(final_contig_data))) {
    if (as.character(final_contig_data$VC[j] == current_cluster)) {
      # print(final_contig_data$Genome[j])
      if (final_contig_data$Dehydration_Status[j] == 1) {
        
        # Append values
        mild_contigs <- c(mild_contigs, final_contig_data$Genome[j])
        mild_ratios <- c(mild_ratios, as.numeric(final_contig_data$`prophage-host_ratio`[j]))
        mild_active <- c(mild_active, as.character(final_contig_data$active[j]))
        
        
        # Append contig (final_contig_data$Genome[j]) to classification_target$mild_contigs[i] array
        # Append prophage-host_ratio (final_contig_data$`prophage-host_ratio`[j]) to classification_target$mild_ratios[i] array
        # Append active status (final_contig_data$active[j]) to classification_target$mild_active[i] array
      } else {
        nonmild_contigs <- c(nonmild_contigs, final_contig_data$Genome[j])
        nonmild_ratios <- c(nonmild_ratios, as.numeric(final_contig_data$`prophage-host_ratio`[j]))
        nonmild_active <- c(nonmild_active, as.character(final_contig_data$active[j]))
        
        # Append contig (final_contig_data$Genome[j]) to classification_target$nonmild_contigs[i] array
        # Append prophage-host_ratio (final_contig_data$`prophage-host_ratio`[j]) to classification_target$nonmild_ratios[i] array
        # Append active status (final_contig_data$active[j]) to classification_target$nonmild_active[i] array
      }
    }
  }
  # After inner loop is done, compute averages for each VC here:
  
  # classification_target$avg_mild_ratios <- the average of the values in mild_ratios
  if (length(mild_ratios) > 0) {
    classification_target$avg_mild_ratios[i] <- mean(mild_ratios, na.rm = TRUE)
  } else {
    classification_target$avg_mild_ratios[i] <- NA
  }
  
  # classification_target$avg_mild_active <- the proportion of the array that is "active" to any other string in mild_active
  if (length(mild_active) > 0) {
    classification_target$avg_mild_active <- sum(mild_active == "active") / length(mild_active)
  } else {
    classification_target$avg_mild_active <- NA
  }
  
  # classification_target$avg_nonmild_ratios <- the average of the values in nonmild_ratios
  if (length(nonmild_ratios) > 0) {
    classification_target$avg_nonmild_ratios[i] <- mean(nonmild_ratios, na.rm = TRUE)
  } else {
    classification_target$avg_nonmild_ratios[i] <- NA
  }
  
  # classification_target$avg_nonmild_active <- the proportion of the array that is "active" to any other string in nonmild_active
  if (length(nonmild_active) > 0) {
    classification_target$avg_nonmild_active <- sum(nonmild_active == "active") / length(nonmild_active)
  } else {
    classification_target$avg_nonmild_active <- NA
  }
}
```

Plot average prophage-host_ratio in mild and non-mild
```{r}
# First filter out rows with NA values in any of the ratio columns
filtered_classification_target <- classification_target %>%
  filter(!is.na(avg_mild_ratios))
filtered_classification_target <- filtered_classification_target %>%
  filter(!is.na(avg_nonmild_ratios))

filtered_classification_target$difference <- 
  (filtered_classification_target$avg_nonmild_ratio - filtered_classification_target$avg_mild_ratio)
filtered_classification_target <- filtered_classification_target %>%
  arrange(desc(difference))

ggplot(filtered_classification_target, aes(x = reorder(VC, -difference), y = difference, fill = difference > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue")) +
  labs(x = "VC", y = "Difference", title = "Mean difference in prophage-host read ratio, NonMild - Mild (per VC)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2))

# Significance testing

median(filtered_classification_target$difference)
# Test normality -> Not normal
shapiro.test(filtered_classification_target$difference)

# Wilcox
wilcox.test(filtered_classification_target$difference, mu = 0, alternative = "two.sided")

# Sign test
pos <- sum(filtered_classification_target$difference > 0)
neg <- sum(filtered_classification_target$difference < 0)
binom.test(pos, pos+neg, p = 0.5)

# Spearman Rank correlation test
# cor.test(filtered_classification_target$difference)
```

Severe vs NonSevere
```{r}
# Make df to store all VCs (that are present in at least two levels of dehydration)
VC_vs_numSamples_tmp <- VC_vs_numSamples %>% select(-`VC Size`)
VC_vs_numSamples_tmp <- VC_vs_numSamples_tmp %>%
  mutate(NonMild = Moderate + Severe)
VC_vs_numSamples_tmp <- VC_vs_numSamples_tmp %>%
  mutate(NonSevere = Moderate + Mild)
classification_target <- merge(VC_vs_numSamples_tmp, stratified_clusters, by = 'VC')

# Then need to initialize columns for the average value
classification_target$avg_severe_ratios    <- rep(NA, nrow(classification_target))
classification_target$avg_severe_active    <- rep(NA, nrow(classification_target))
classification_target$avg_nonsevere_ratios <- rep(NA, nrow(classification_target))
classification_target$avg_nonsevere_active <- rep(NA, nrow(classification_target))

for (i in seq_len(nrow(classification_target))) {
  current_cluster <- as.character(classification_target$VC[i])
  
  severe_contigs <- character(0)
  severe_ratios <- numeric(0)
  severe_active <- character(0)
  
  nonsevere_contigs <- character(0)
  nonsevere_ratios <- numeric(0)
  nonsevere_active <- character(0)

  for (j in seq_len(nrow(final_contig_data))) {
    if (as.character(final_contig_data$VC[j] == current_cluster)) {
      if (final_contig_data$Dehydration_Status[j] == 3) {
        
        # Append values
        severe_contigs <- c(severe_contigs, final_contig_data$Genome[j])
        severe_ratios <- c(severe_ratios, as.numeric(final_contig_data$`prophage-host_ratio`[j]))
        severe_active <- c(severe_active, as.character(final_contig_data$active[j]))
        
      } else {
        nonsevere_contigs <- c(nonsevere_contigs, final_contig_data$Genome[j])
        nonsevere_ratios <- c(nonsevere_ratios, as.numeric(final_contig_data$`prophage-host_ratio`[j]))
        nonsevere_active <- c(nonsevere_active, as.character(final_contig_data$active[j]))
        
      }
    }
  }
  # After inner loop is done, compute averages for each VC here:
  
  # classification_target$avg_mild_ratios <- the average of the values in mild_ratios
  if (length(severe_ratios) > 0) {
    classification_target$avg_severe_ratios[i] <- mean(severe_ratios, na.rm = TRUE)
  } else {
    classification_target$avg_severe_ratios[i] <- NA
  }
  
  # classification_target$avg_mild_active <- the proportion of the array that is "active" to any other string in mild_active
  if (length(severe_active) > 0) {
    classification_target$avg_severe_active <- sum(severe_active == "active") / length(severe_active)
  } else {
    classification_target$avg_severe_active <- NA
  }
  
  # classification_target$avg_nonmild_ratios <- the average of the values in nonmild_ratios
  if (length(nonsevere_ratios) > 0) {
    classification_target$avg_nonsevere_ratios[i] <- mean(nonsevere_ratios, na.rm = TRUE)
  } else {
    classification_target$avg_nonsevere_ratios[i] <- NA
  }
  
  # classification_target$avg_nonmild_active <- the proportion of the array that is "active" to any other string in nonmild_active
  if (length(nonsevere_active) > 0) {
    classification_target$avg_nonsevere_active <- sum(nonsevere_active == "active") / length(nonsevere_active)
  } else {
    classification_target$avg_nonsevere_active <- NA
  }
}
```

```{r}
# First filter out rows with NA values in any of the ratio columns
filtered_classification_target <- classification_target %>%
  filter(!is.na(avg_severe_ratios))
filtered_classification_target <- filtered_classification_target %>%
  filter(!is.na(avg_nonsevere_ratios))

filtered_classification_target$difference <- 
  (filtered_classification_target$avg_severe_ratio - filtered_classification_target$avg_nonsevere_ratio)
filtered_classification_target <- filtered_classification_target %>%
  arrange(desc(difference))

ggplot(filtered_classification_target, aes(x = reorder(VC, -difference), y = difference, fill = difference > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue")) +
  labs(x = "VC", y = "Difference", title = "Mean difference in prophage-host read ratio, severe - nonSevere (per VC)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2))

# Significance testing

median(filtered_classification_target$difference)

# Test normality -> Not normal
shapiro.test(filtered_classification_target$difference)

# Wilcox -> p = 0.003635 
wilcox.test(filtered_classification_target$difference, mu = 0, alternative = "two.sided")

```

Within these VCs what are the taxa?
```{r}

```


## Antibiotics
```{r}
# Add antibiotic data associated with each contig's sample
for (i in seq_len(nrow(final_contig_data))) {
  for (j in seq_len(nrow(MG_data_NM))) {
    if (as.character(MG_data_NM$Sample[j]) == as.character(final_contig_data$Sample[i])) {
      final_contig_data$DOX[i] <- MG_data_NM$DOX[j]
      final_contig_data$CIP[i] <- MG_data_NM$CIP[j]
      final_contig_data$AZI[i] <- MG_data_NM$AZI[j]
    } 
  }
}
```

```{r}
abx_VC_vs_numSamples <- final_contig_data %>%
  group_by(VC) %>%
  summarise(
    DOX_pos = sum(DOX >= 0.13, na.rm = TRUE),
    DOX_neg = sum(DOX < 0.13, na.rm = TRUE),
    CIP_pos = sum(CIP >= 0.063, na.rm = TRUE),
    CIP_neg = sum(CIP < 0.063, na.rm = TRUE),
    AZI_pos = sum(AZI >= 8, na.rm = TRUE),
    AZI_neg = sum(AZI < 8, na.rm = TRUE)
  ) %>% 
  ungroup()

# Keep only row where stratification in DOX exposure exists
print(nrow(abx_VC_vs_numSamples))
DOX_VC_vs_numSamples <- abx_VC_vs_numSamples %>%
  filter((!abx_VC_vs_numSamples$DOX_pos == 0 & !abx_VC_vs_numSamples$DOX_neg == 0))
print(nrow(DOX_VC_vs_numSamples))
# 65 remaining clusters

# Keep only row where stratification in CIP exposure exists
print(nrow(abx_VC_vs_numSamples))
CIP_VC_vs_numSamples <- abx_VC_vs_numSamples %>%
  filter((!abx_VC_vs_numSamples$CIP_pos == 0 & !abx_VC_vs_numSamples$CIP_neg == 0))
print(nrow(CIP_VC_vs_numSamples))
# 288 remaining clusters

# Keep only row where stratification in AZI exposure exists
print(nrow(abx_VC_vs_numSamples))
AZI_VC_vs_numSamples <- abx_VC_vs_numSamples %>%
  filter((!abx_VC_vs_numSamples$AZI_pos == 0 & !abx_VC_vs_numSamples$AZI_neg == 0))
print(nrow(AZI_VC_vs_numSamples))
# 134 remaining clusters
```

```{r}
#abx_classification_target <- AZI_VC_vs_numSamples
#antibiotic_target <- "AZI"
#threshold <- 8

#abx_classification_target <- CIP_VC_vs_numSamples
#antibiotic_target <- "CIP"
#threshold <- 0.063

abx_classification_target <- DOX_VC_vs_numSamples
antibiotic_target <- "DOX"
threshold <- 0.13

# Then need to initialize columns for the average value
abx_classification_target$avg_pos_ratios <- rep(NA, nrow(abx_classification_target))
abx_classification_target$avg_pos_active <- rep(NA, nrow(abx_classification_target))
abx_classification_target$avg_neg_ratios <- rep(NA, nrow(abx_classification_target))
abx_classification_target$avg_neg_active <- rep(NA, nrow(abx_classification_target))

for (i in seq_len(nrow(abx_classification_target))) {
  current_cluster <- as.character(abx_classification_target$VC[i])
  
  pos_ratios <- numeric(0)
  pos_active <- character(0)
  
  neg_ratios <- numeric(0)
  neg_active <- character(0)

  for (j in seq_len(nrow(final_contig_data))) {
    if (as.character(final_contig_data$VC[j] == current_cluster)) {
      # print(final_contig_data$Genome[j])
      if (!is.na(final_contig_data[[antibiotic_target]][j]) && final_contig_data[[antibiotic_target]][j] >= threshold) {
        
        # Append values
        pos_ratios <- c(pos_ratios, as.numeric(final_contig_data$`prophage-host_ratio`[j]))
        pos_active <- c(pos_active, as.character(final_contig_data$active[j]))
        
      } else {
        neg_ratios <- c(neg_ratios, as.numeric(final_contig_data$`prophage-host_ratio`[j]))
        neg_active <- c(neg_active, as.character(final_contig_data$active[j]))
      }
    }
  }
  # After inner loop is done, compute averages for each VC here:
  
  # classification_target$avg_mild_ratios <- the average of the values in mild_ratios
  if (length(pos_ratios) > 0) {
    abx_classification_target$avg_pos_ratios[i] <- mean(pos_ratios, na.rm = TRUE)
  } else {
    abx_classification_target$avg_pos_ratios[i] <- NA
  }
  
  # classification_target$avg_mild_active <- the proportion of the array that is "active" to any other string in mild_active
  if (length(pos_active) > 0) {
    abx_classification_target$avg_pos_active <- sum(pos_active == "active") / length(pos_active)
  } else {
    abx_classification_target$avg_pos_active <- NA
  }
  
  # classification_target$avg_nonmild_ratios <- the average of the values in nonmild_ratios
  if (length(neg_ratios) > 0) {
    abx_classification_target$avg_neg_ratios[i] <- mean(neg_ratios, na.rm = TRUE)
  } else {
    abx_classification_target$avg_neg_ratios[i] <- NA
  }
  
  # classification_target$avg_nonmild_active <- the proportion of the array that is "active" to any other string in nonmild_active
  if (length(neg_active) > 0) {
    abx_classification_target$avg_neg_active <- sum(neg_active == "active") / length(neg_active)
  } else {
    abx_classification_target$avg_neg_active <- NA
  }
}

# Filter out rows with NA values
abx_classification_target <- abx_classification_target %>%
  filter(!is.na(avg_pos_ratios))
abx_classification_target <- abx_classification_target %>%
  filter(!is.na(avg_neg_ratios))

# Compute difference 
abx_classification_target$difference <-
  (abx_classification_target$avg_pos_ratios - abx_classification_target$avg_neg_ratios)
abx_classification_target <- abx_classification_target %>%
  arrange(desc(difference))

# Plot
ggplot(abx_classification_target, aes(x = reorder(VC, -difference), y = difference, fill = difference > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue")) +
  labs(x = "VC", y = "Difference", title = "(DOX) Drug Presence - Drug Negative (per VC)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2))

# Significance testing

median(abx_classification_target$difference)

# Test normality -> Not normal
shapiro.test(abx_classification_target$difference)

# Wilcox -> p = 0.003635 
wilcox.test(abx_classification_target$difference, mu = 0, alternative = "two.sided")
```

## Separate data by dehydration status
It is possible that dehydration status and antibiotic exposure are collinear and 
"cancelling each other out". To control for collinearity between both variables,
I will investigate the effect of prophage induction 

```{r}
# Subset the data into only severe
filtered_genome_by_genome <- genome_by_genome_overview %>% 
  filter(grepl("^k", Genome))
clustered_genome_by_genome <- filtered_genome_by_genome %>% 
  filter(`VC Status` == "Clustered")

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
```

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

VC_vs_numSamples <- merge(VC_vs_numSamples, vclusters, by = 'VC')

# unstratified_clusters <- sum(rowSums(VC_vs_numSamples[, c("Mild", "Moderate", "Severe")] > 0) == 1)
# print(unstratified_clusters)

propagate_by_contigs <- merge(sample_by_contig, propagate_results_by_prophage_liberal, by = 'Genome')

# Merge sample dehydration status to filtered_genome_by_genome
sample_and_status <- X00_metadata %>%
  select(Sample, Dehydration_Status)

# Merge filtered genome with sample_and_status, then remove duplicate columns, active and prophage-host_ratio


filtered_genome_by_genome <- merge(filtered_genome_by_genome, propagate_by_contigs, by = 'Genome')
filtered_genome_by_genome <- merge(filtered_genome_by_genome, sample_and_status, by = 'Sample')
filtered_genome_by_genome <- filtered_genome_by_genome %>% select(-active, -`prophage-host_ratio`)
filtered_genome_by_genome <- filtered_genome_by_genome %>%
  filter(Dehydration_Status == 3)


final_contig_data <- merge(filtered_genome_by_genome, propagate_results_by_prophage_liberal, by = 'Genome')
final_contig_data <- final_contig_data %>%
  select(
    Genome, 
    VC,
    Sample, 
    Dehydration_Status,
    active, 
    `prophage-host_ratio`
  )
final_contig_data <- final_contig_data %>%
  filter(!is.na(VC))

# Add antibiotic data associated with each contig's sample
for (i in seq_len(nrow(final_contig_data))) {
  for (j in seq_len(nrow(MG_data_NM))) {
    if (as.character(MG_data_NM$Sample[j]) == as.character(final_contig_data$Sample[i])) {
      final_contig_data$DOX[i] <- MG_data_NM$DOX[j]
      final_contig_data$CIP[i] <- MG_data_NM$CIP[j]
      final_contig_data$AZI[i] <- MG_data_NM$AZI[j]
    } 
  }
}
```

```{r}
abx_VC_vs_numSamples <- final_contig_data %>%
  group_by(VC) %>%
  summarise(
    DOX_pos = sum(DOX >= 0.13, na.rm = TRUE),
    DOX_neg = sum(DOX < 0.13, na.rm = TRUE),
    CIP_pos = sum(CIP >= 0.063, na.rm = TRUE),
    CIP_neg = sum(CIP < 0.063, na.rm = TRUE),
    AZI_pos = sum(AZI >= 8, na.rm = TRUE),
    AZI_neg = sum(AZI < 8, na.rm = TRUE)
  ) %>% 
  ungroup()

# Keep only row where stratification in DOX exposure exists
print(nrow(abx_VC_vs_numSamples))
DOX_VC_vs_numSamples <- abx_VC_vs_numSamples %>%
  filter((!abx_VC_vs_numSamples$DOX_pos == 0 & !abx_VC_vs_numSamples$DOX_neg == 0))
print(nrow(DOX_VC_vs_numSamples))
# 65 remaining clusters

# Keep only row where stratification in CIP exposure exists
print(nrow(abx_VC_vs_numSamples))
CIP_VC_vs_numSamples <- abx_VC_vs_numSamples %>%
  filter((!abx_VC_vs_numSamples$CIP_pos == 0 & !abx_VC_vs_numSamples$CIP_neg == 0))
print(nrow(CIP_VC_vs_numSamples))
# 288 remaining clusters

# Keep only row where stratification in AZI exposure exists
print(nrow(abx_VC_vs_numSamples))
AZI_VC_vs_numSamples <- abx_VC_vs_numSamples %>%
  filter((!abx_VC_vs_numSamples$AZI_pos == 0 & !abx_VC_vs_numSamples$AZI_neg == 0))
print(nrow(AZI_VC_vs_numSamples))
# 134 remaining clusters
```

```{r}
abx_classification_target <- AZI_VC_vs_numSamples
antibiotic_target <- "AZI"
threshold <- 8

#abx_classification_target <- CIP_VC_vs_numSamples
#antibiotic_target <- "CIP"
#threshold <- 0.063

#abx_classification_target <- DOX_VC_vs_numSamples
#antibiotic_target <- "DOX"
#threshold <- 0.13

# Then need to initialize columns for the average value
abx_classification_target$avg_pos_ratios <- rep(NA, nrow(abx_classification_target))
abx_classification_target$avg_pos_active <- rep(NA, nrow(abx_classification_target))
abx_classification_target$avg_neg_ratios <- rep(NA, nrow(abx_classification_target))
abx_classification_target$avg_neg_active <- rep(NA, nrow(abx_classification_target))

for (i in seq_len(nrow(abx_classification_target))) {
  current_cluster <- as.character(abx_classification_target$VC[i])
  
  pos_ratios <- numeric(0)
  pos_active <- character(0)
  
  neg_ratios <- numeric(0)
  neg_active <- character(0)

  for (j in seq_len(nrow(final_contig_data))) {
    if (as.character(final_contig_data$VC[j] == current_cluster)) {
      # print(final_contig_data$Genome[j])
      if (!is.na(final_contig_data[[antibiotic_target]][j]) && final_contig_data[[antibiotic_target]][j] >= threshold) {
        
        # Append values
        pos_ratios <- c(pos_ratios, as.numeric(final_contig_data$`prophage-host_ratio`[j]))
        pos_active <- c(pos_active, as.character(final_contig_data$active[j]))
        
      } else {
        neg_ratios <- c(neg_ratios, as.numeric(final_contig_data$`prophage-host_ratio`[j]))
        neg_active <- c(neg_active, as.character(final_contig_data$active[j]))
      }
    }
  }
  # After inner loop is done, compute averages for each VC here:
  
  # classification_target$avg_mild_ratios <- the average of the values in mild_ratios
  if (length(pos_ratios) > 0) {
    abx_classification_target$avg_pos_ratios[i] <- mean(pos_ratios, na.rm = TRUE)
  } else {
    abx_classification_target$avg_pos_ratios[i] <- NA
  }
  
  # classification_target$avg_mild_active <- the proportion of the array that is "active" to any other string in mild_active
  if (length(pos_active) > 0) {
    abx_classification_target$avg_pos_active <- sum(pos_active == "active") / length(pos_active)
  } else {
    abx_classification_target$avg_pos_active <- NA
  }
  
  # classification_target$avg_nonmild_ratios <- the average of the values in nonmild_ratios
  if (length(neg_ratios) > 0) {
    abx_classification_target$avg_neg_ratios[i] <- mean(neg_ratios, na.rm = TRUE)
  } else {
    abx_classification_target$avg_neg_ratios[i] <- NA
  }
  
  # classification_target$avg_nonmild_active <- the proportion of the array that is "active" to any other string in nonmild_active
  if (length(neg_active) > 0) {
    abx_classification_target$avg_neg_active <- sum(neg_active == "active") / length(neg_active)
  } else {
    abx_classification_target$avg_neg_active <- NA
  }
}

# Filter out rows with NA values
abx_classification_target <- abx_classification_target %>%
  filter(!is.na(avg_pos_ratios))
abx_classification_target <- abx_classification_target %>%
  filter(!is.na(avg_neg_ratios))

# Compute difference 
abx_classification_target$difference <-
  (abx_classification_target$avg_pos_ratios - abx_classification_target$avg_neg_ratios)
abx_classification_target <- abx_classification_target %>%
  arrange(desc(difference))

# Plot
ggplot(abx_classification_target, aes(x = reorder(VC, -difference), y = difference, fill = difference > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue")) +
  labs(x = "VC", y = "Difference", title = "(AZI) Drug Presence - Drug Negative (per VC)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2))

# Significance testing

median(abx_classification_target$difference)

# Test normality -> Not normal
shapiro.test(abx_classification_target$difference)

# Wilcox -> p = 0.003635 
wilcox.test(abx_classification_target$difference, mu = 0, alternative = "two.sided")
```

## GLMM

### Libraries
```{r}
library(Matrix)
require(lme4)
#library(ordinal)
library(performance)
```

```{r}
glmm_propagate_results <- propagate_by_contigs %>%
  filter(!is.na(`prophage-host_ratio`))

for (sample in seq_len(nrow(X00_metadata))) {
  current_ratios <- numeric(0)
  for (contig in seq_len(nrow(glmm_propagate_results))) {
    if (as.character(glmm_propagate_results$Sample[contig]) == as.character(X00_metadata$Sample[sample])) {
      current_ratios <- c(current_ratios, as.numeric(glmm_propagate_results$`prophage-host_ratio`[contig]))
    }
  }
  result <- mean(current_ratios)
  X00_metadata$`prophage-host_ratio`[sample] <- result
} 

glmm_data <- merge(X00_metadata, MG_data_NM, by = "Sample")
glmm_data <- glmm_data %>%
  select(Sample,
         Dehydration_Status,
         `prophage-host_ratio`,
         Area_Code,
         AZI,
         CIP, 
         DOX)
glmm_data$Dehydration_Status <- factor(glmm_data$Dehydration_Status, 
                                       ordered = TRUE)
```

```{r}
# glmm_model <- clmm(Dehydration_Status ~ `prophage-host_ratio` + CIP + (1 | Area_Code), data = glmm_data)
glmm_model <- glm(Dehydration_Status~`prophage-host_ratio` + CIP + (1 | Area_Code), data = glmm_data, family = poisson)
#check_collinearity(glmm_model)
```
