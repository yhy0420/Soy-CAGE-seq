 Soy-CAGE-seq
## Scripts for analyzing soybean nanoCAGE data

'01.CAGE-alignment.sh' is used to process nanoCAGE-seq data, including adapter removal, rRNA, chloroplast and mitochondrial read filtering, and reads mapping. The final output is a BAM file, which is used for subsequent TSS (transcription start site) identification.

'02.CAGEr.rmd' uses the CAGEr package to identify transcription start sites (TSSs) and profiling.

'03.anno_Soy_TSS.R' complements the published genome annotation with TSSs identified from nanoCAGE-seq data.

'04.shape.py' is used to process promoter region and identified TSSs, and classify promoters based on information entropy (Shannon index, SI).

The four bw files are normalized reads depth profiling in bigwig format, which could be virtualized using tools such as [IGV](https://igv.org/). 

Root_TC.bed and Shoot_TC.bed are the annotated CTSSs.

Root_TC.gff3 and Shoot_TC.gff3 are the annotated dominant CTSSs.
