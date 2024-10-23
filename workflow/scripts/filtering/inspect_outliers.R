# Setup
suppressMessages(library(tidyverse))

# Input
nucc <- read.csv(file=snakemake@input[["nucc"]], header=TRUE, check.names=FALSE,
                 stringsAsFactors=FALSE, na.strings=c('', 'NA', -1, '-1')) %>%
  mutate(NUCCORE_Topology = ifelse(is.na(NUCCORE_Topology), 'not_set', NUCCORE_Topology ))
typing <- read.csv(file=snakemake@input[["typing"]], header=TRUE, check.names=FALSE,
                 stringsAsFactors=FALSE, na.strings=c('', 'NA', -1, '-1', '-'))
amr <- read.csv(file=snakemake@input[["amr"]], header=TRUE, check.names=FALSE, sep="\t",
                 stringsAsFactors=FALSE, na.strings=c('', 'NA', -1, '-1'))
prot <- read.csv(file=snakemake@input[["proteins"]], header=TRUE, check.names=FALSE,
                 stringsAsFactors=FALSE, na.strings=c('', 'NA', -1, '-1'))

# Inspect GC and Length outliers
gc_outliers <- nucc %>% filter(
  NUCCORE_GC < quantile(nucc$NUCCORE_GC, 0.25)-1.5*IQR(nucc$NUCCORE_GC) | NUCCORE_GC > quantile(nucc$NUCCORE_GC, 0.75)+1.5*IQR(nucc$NUCCORE_GC)) %>% pull(NUCCORE_GC)
len_outliers <- nucc %>% filter(
  log10(NUCCORE_Length) < quantile(log10(nucc$NUCCORE_Length), 0.25)-1.5*IQR(log10(nucc$NUCCORE_Length)) | log10(NUCCORE_Length) > quantile(log10(NUCCORE_Length), 0.75)+1.5*IQR(log10(NUCCORE_Length))) %>% pull(NUCCORE_Length)


to_eval <- nucc %>%
  filter(NUCCORE_GC %in% gc_outliers | NUCCORE_Length %in% len_outliers) %>%
  left_join(typing %>% select(NUCCORE_ACC, 
                              `rep_type(s)`, `relaxase_type(s)`, mpf_type,
                              `orit_type(s)`, PMLST_scheme, PMLST_sequence_type), 
            by=join_by(NUCCORE_ACC)) %>%
  left_join(amr %>% group_by(NUCCORE_ACC) %>% count(name="AMR_genes_count"), by=join_by(NUCCORE_ACC)) %>%
  left_join(prot %>% group_by(NUCCORE_ACC) %>% count(name="elements_count"), by=join_by(NUCCORE_ACC)) %>%
  select(NUCCORE_ACC, NUCCORE_Length, NUCCORE_Completeness,
         NUCCORE_GC,
         `rep_type(s)`, `relaxase_type(s)`, mpf_type,
         `orit_type(s)`, PMLST_scheme, PMLST_sequence_type,
         VIRAL_category, VIRAL_putative,
         AMR_genes_count, elements_count
         )
print(paste("# Candidates to evalutate:", nrow(to_eval)))

# Discard if any typing present, neither AMR genes or annotations
to_discard <- to_eval %>%
  filter_at(vars(
    `rep_type(s)`, `relaxase_type(s)`, mpf_type,
    `orit_type(s)`, PMLST_scheme, PMLST_sequence_type,
    AMR_genes_count, elements_count
    ), all_vars(is.na(.))) %>%
  mutate(STATUS = "suppressed") %>%
  mutate(STATUS_note = "Anomalous GC content/Length: Absence of any typing (replicon, relaxases, MPF, oriT, pMLST) and genetic elements (i.e. protein-coding genes)") %>% 
  select(NUCCORE_ACC, STATUS, STATUS_note)

print(paste("# Candidates to discard:", nrow(to_eval)))

length_cutoff <- nucc %>%
  filter(!(NUCCORE_ACC %in% to_discard$NUCCORE_ACC)) %>%
  filter(NUCCORE_Length > 4000000 | NUCCORE_Length < 1000) %>%
  mutate(STATUS = "suppressed") %>%
  mutate(STATUS_note = "Anomalous GC content/Length: Length > 4MB or Length < 1Kbp") %>% 
  select(NUCCORE_ACC, STATUS, STATUS_note)

print(paste("# Candidates discarded from length cutoff:", nrow(length_cutoff)))
# Inspect others
to_discard <- rbind(to_discard, length_cutoff) 


to_discard <- rbind(to_discard, 
data.frame(
  NUCCORE_ACC=c("NZ_CP126587.1"),
  STATUS=c("suppressed"),
  STATUS_note=c("Anomalous GC content/Length: artifact region")
  ) )

print(paste("# Plasmids to discard", nrow(to_discard)))
print(table(to_discard$STATUS_note))


write.csv(to_discard, snakemake@output[['csv']], row.names=FALSE)



