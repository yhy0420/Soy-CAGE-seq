setwd("/home/dell/file/SangLab/CAGE-seq/250528/TC7/GFF/")
rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

# 读取并处理 Shoot.gff3 文件
shoot_S <- read_tsv('Shoot_TC7.gff3', comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = sub(";tpm=.*", "", annotation),  # 提取基因ID
    tpm = as.numeric(sub(".*;tpm=", "", annotation))  # 提取TPM值并转换为数值型
  )

# 读取并处理 Root.gff3 文件
root_S <- read_tsv('Root_TC7.gff3', comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = sub(";tpm=.*", "", annotation),  # 提取基因ID
    tpm = as.numeric(sub(".*;tpm=", "", annotation))  # 提取TPM值并转换为数值型
  )

# 读取并处理 Glymax_cage.gff 文件
reference_S <- read_tsv('../../Glymax_cage.gff', comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = substr(annotation, 4, 18)
  )

# 保留TPM值较大的行
Shoot_S <- shoot_S %>%
  group_by(gene) %>%
  slice_max(order_by = tpm, n = 1) %>%
  ungroup()  

Root_S <- root_S %>%
  group_by(gene) %>%
  slice_max(order_by = tpm, n = 1) %>%
  ungroup()  

# 合并 shoot 和 root dTSSs
correctedTSS_S <- union(Shoot_S, Root_S)

# 保留TPM值较大的行
correctedTSS_S <- correctedTSS_S %>%
  group_by(gene) %>%
  slice_max(order_by = tpm, n = 1) %>%
  ungroup()  

correctedTSS_S <- correctedTSS_S %>% select(-tpm)

# 获取参考 TSSs
referenceTSS_S <- reference_S %>%
  mutate(
    TSS = if_else(strand == '+', start, stop)
  ) %>%
  select(gene, TSS) %>%
  deframe()

# 生成 multipleTSS
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

# 生成 correctedTSS1
correctedTSS1 <- correctedTSS_S %>%
  mutate(
    gene = if_else(annotation %in% names(multipleTSS), multipleTSS[annotation], gene)
  ) %>%
  select(gene, stop) %>%
  deframe()
# export names of genes with a dTSS in the Grotewold data
write_lines(names(correctedTSS1), 'Soy_corrected_TSS_TC7.txt')

# 修正注释时只替换参考基因组上最接近的位点
TSS_S1 <- reference_S %>%
  mutate(
    is_gene_or_five_prime_utr = type %in% c('gene', 'five_prime_UTR'),
    start = if_else(strand == '+' & gene %in% names(correctedTSS1) & is_gene_or_five_prime_utr, correctedTSS1[gene], start),
    stop = if_else(strand == '-' & gene %in% names(correctedTSS1) & is_gene_or_five_prime_utr, correctedTSS1[gene], stop)
  ) %>% 
  select(-gene, -is_gene_or_five_prime_utr)
write_tsv(TSS_S1, 'Soy_protein_coding_genes_TC7_1.gff3', col_names = FALSE)

# 修正 TSS 注释 gene, five_prime_UTR, transcript
TSS_S2 <- reference_S %>%
  mutate(
    is_gene_or_five_prime_utr_or_transcript = type %in% c('gene', 'five_prime_UTR', 'transcript'),
    start = if_else(strand == '+' & gene %in% names(correctedTSS1) & is_gene_or_five_prime_utr_or_transcript, correctedTSS1[gene], start),
    stop = if_else(strand == '-' & gene %in% names(correctedTSS1) & is_gene_or_five_prime_utr_or_transcript, correctedTSS1[gene], stop)
  ) %>% 
  select(-gene, -is_gene_or_five_prime_utr_or_transcript)


# 保存修正后的注释到文件
write_tsv(TSS_S2, 'Soy_protein_coding_genes_TC7_2.gff3', col_names = FALSE)
