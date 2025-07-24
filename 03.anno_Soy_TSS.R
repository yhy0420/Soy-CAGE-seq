# CAGE anno
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

shoot_S <- read_tsv('Shoot_TC7.gff3', comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = sub(";tpm=.*", "", annotation), 
    tpm = as.numeric(sub(".*;tpm=", "", annotation))
  )

root_S <- read_tsv('Root_TC7.gff3', comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = sub(";tpm=.*", "", annotation), 
    tpm = as.numeric(sub(".*;tpm=", "", annotation))
  )

reference_S <- read_tsv('../../Glymax_cage.gff', comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = substr(annotation, 4, 18)
  )

Shoot_S <- shoot_S %>%
  group_by(gene) %>%
  slice_max(order_by = tpm, n = 1) %>%
  ungroup()  

Root_S <- root_S %>%
  group_by(gene) %>%
  slice_max(order_by = tpm, n = 1) %>%
  ungroup()  

correctedTSS_S <- union(Shoot_S, Root_S)

correctedTSS_S <- correctedTSS_S %>%
  group_by(gene) %>%
  slice_max(order_by = tpm, n = 1) %>%
  ungroup()  

correctedTSS_S <- correctedTSS_S %>% select(-tpm)

referenceTSS_S <- reference_S %>%
  mutate(
    TSS = if_else(strand == '+', start, stop)
  ) %>%
  select(gene, TSS) %>%
  deframe()

multipleTSS <- correctedTSS_S %>%
  mutate(
    refTSS = referenceTSS_S[gene]
  ) %>%
  group_by(across(c(-gene, -refTSS))) %>%
  summarise(
    gene = first(gene[abs(stop - refTSS) == min(abs(stop - refTSS))])
  ) %>%
  ungroup() %>%
  select(annotation, gene) %>%
  deframe()

correctedTSS1 <- correctedTSS_S %>%
  mutate(
    gene = if_else(annotation %in% names(multipleTSS), multipleTSS[annotation], gene)
  ) %>%
  select(gene, stop) %>%
  deframe()
# export names of genes with a dTSS in the Grotewold data
write_lines(names(correctedTSS1), 'Soy_corrected_TSS_TC.txt')

TSS_S1 <- reference_S %>%
  mutate(
    is_gene_or_five_prime_utr = type %in% c('gene', 'five_prime_UTR'),
    start = if_else(strand == '+' & gene %in% names(correctedTSS1) & is_gene_or_five_prime_utr, correctedTSS1[gene], start),
    stop = if_else(strand == '-' & gene %in% names(correctedTSS1) & is_gene_or_five_prime_utr, correctedTSS1[gene], stop)
  ) %>% 
  select(-gene, -is_gene_or_five_prime_utr)
write_tsv(TSS_S1, 'Soy_protein_coding_genes_TC.gff3', col_names = FALSE)
