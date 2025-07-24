#!/bin/bash

# Input/Output settings
SAMPLE_LIST="../sample.list"
RAW_DATA_DIR="/path/to/raw_fastq"
OUTPUT_DIR="$(pwd)/CAGE_processed"
mkdir -p ${OUTPUT_DIR}/{1.raw_extract,2.adapter_trim,3.rRNA_filter,4.genome_align,5.final_bams}

# Sequence definitions
ADAPTER_1="TATAGGG"
ADAPTER_2="AGATCGGAAGAGC" 
PRIMER_SEQ="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT" # Primer sequence

# Reference genomes
RRNA_DB="/path/to/rRNA_index"
GENOME_REF="/path/to/genome_index"
CHR_TO_REMOVE="NC_007942.1,NC_020455.1" # Chloroplast and mitochondria IDs

# ==================== PROCESSING FUNCTION ====================
process_sample() {
    local sample=$1
    local log_file="${OUTPUT_DIR}/${sample}.log"
    
    echo "[$(date)] Processing ${sample}" > ${log_file}
    local start_time=$(date +%s)

    # ----- STEP 1: Extract reads with both adapters -----
    echo "[STEP1] Extracting reads with adapters" >> ${log_file}
    local raw_input="${RAW_DATA_DIR}/${sample}_R1.fastq.gz"
    local step1_out="${OUTPUT_DIR}/1.raw_extract/${sample}_with_adapters.fastq.gz"
    
    zcat ${raw_input} | awk -v ad5="${ADAPTER_1}" -v ad3="${ADAPTER_2}" '
        BEGIN { RS = "\n@"; FS = "\n" }
        NR > 1 && $2 ~ ad5 && $2 ~ ad3 {
            print "@" $1 "\n" $2 "\n" $3 "\n" $4
        }
    ' | pigz -p 10 > ${step1_out} 2>> ${log_file}

    if [[ ! -s ${step1_out} ]]; then
        echo "[ERROR] No adapter-containing reads found!" >> ${log_file}
        return 1
    fi

    # ----- STEP 2: Adapter/primer trimming -----
    echo "[STEP2] Trimming adapters/primers" >> ${log_file}
    local step2_out="${OUTPUT_DIR}/2.adapter_trim/${sample}_clean.fastq.gz"
    
    # Generate primer k-mers
    local primer_kmer=$(for ((i=0; i<=${#PRIMER_SEQ}-10; i++)); do
        echo "${PRIMER_SEQ:i:10}"; 
    done | paste -sd "|")

    cutadapt -a "${ADAPTER_1}...${ADAPTER_2}" \
             -e 0.1 -q 20 -m 15 -j 20 \
             --max-n=0 --max-expected-errors=0.3 \
             -o ${step2_out} ${step1_out} >> ${log_file} 2>&1

    # ----- STEP 3: rRNA removal -----
    echo "[STEP3] Removing rRNA reads" >> ${log_file}
    local step3_out="${OUTPUT_DIR}/3.rRNA_filter/${sample}_non_rrna.fastq.gz"
    
    bowtie2 --local --sensitive-local -x ${RRNA_DB} \
            -U ${step2_out} --un-gz ${step3_out} \
            -S ${OUTPUT_DIR}/3.rRNA_filter/${sample}_rrna.sam \
            --threads 20 --no-unal >> ${log_file} 2>&1

    # ----- STEP 4: Genome alignment -----
    echo "[STEP4] Genome alignment" >> ${log_file}
    local step4_dir="${OUTPUT_DIR}/4.genome_align"
    
    bowtie2 --local --sensitive-local -x ${GENOME_REF} \
            -U ${step3_out} -p 20 \
            -S ${step4_dir}/${sample}.sam >> ${log_file} 2>&1
    
    # Convert and sort BAM
    samtools view -@ 10 -Shu ${step4_dir}/${sample}.sam | \
    samtools sort -@ 10 -o ${step4_dir}/${sample}.sorted.bam -
    samtools index -@ 10 ${step4_dir}/${sample}.sorted.bam

    # ----- STEP 5: Final filtering -----
    echo "[STEP5] Final filtering" >> ${log_file}
    local step5_dir="${OUTPUT_DIR}/5.final_bams"
    
    # Extract high-quality reads
    samtools view -@ 10 -F 4 -q 20 -u ${step4_dir}/${sample}.sorted.bam | \
    samtools sort -@ 10 -o ${step5_dir}/${sample}.highQual.bam -
    
    # Remove organelle reads
    samtools view -@ 10 ${step5_dir}/${sample}.highQual.bam | \
    awk -v ids="${CHR_TO_REMOVE}" 'BEGIN{split(ids,a,","); for(i in a) chr[a[i]]=1} 
                                   $3 in chr {print $1}' | sort -u > ${step5_dir}/${sample}.remove.list
    
    samtools view -@ 10 -h -U ${step5_dir}/${sample}.final.bam \
              -N ${step5_dir}/${sample}.remove.list \
              ${step5_dir}/${sample}.highQual.bam
    
    # Final indexing and stats
    samtools index -@ 10 ${step5_dir}/${sample}.final.bam
    samtools flagstat -@ 10 ${step5_dir}/${sample}.final.bam > ${step5_dir}/${sample}.stats

    # ----- Cleanup and reporting -----
    local end_time=$(date +%s)
    echo "[STATS] Processing completed in $((end_time - start_time)) seconds" >> ${log_file}
    echo "[STATS] Final reads: $(samtools view -c -F 4 ${step5_dir}/${sample}.final.bam)" >> ${log_file}
}

# ==================== MAIN EXECUTION ====================
export -f process_sample
export OUTPUT_DIR RAW_DATA_DIR ADAPTER_5 ADAPTER_3 PRIMER_SEQ
export RRNA_DB GENOME_REF CHR_TO_REMOVE

# Process samples in parallel
parallel -j $SLURM_NTASKS --eta --joblog ${OUTPUT_DIR}/parallel.log \
         "process_sample {}" :::: ${SAMPLE_LIST}

echo "[$(date)] Pipeline completed" | tee -a ${OUTPUT_DIR}/master.log
