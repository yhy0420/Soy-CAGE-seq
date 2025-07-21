# Soy-CAGE-seq
## Scripts for analyzing soybean CAGE data
Data processing and CAGEr analysis

‘CAGE-alignment.sh’ is used to process CAGE-seq data, including adapter removal, rRNA filtering, chloroplast and mitochondrial contamination removal, and genome alignment. The final output is a BAM file, which is used for subsequent TSS (transcription start site) analysis.

‘CAGEr.rmd’ uses the CAGEr package to process CAGE-seq data.

‘anno_Soy_TSS.R’ corrects GFF files based on transcription start sites (TSS) identified from CAGE-seq data.
