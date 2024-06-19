#!/usr/bin/bash

cd /master/nplatt/sch_hae_scan

QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -p -1023" 

mkdir -p results/haplotype_caller/contig_lists

cat data/schHae2_wmito.fa.fai \
    | cut -f1 \
    | shuf \
    | split \
        --additional-suffix .list \
        --numeric-suffixes=1 \
        -l 100 \
        -a 1 \
        - \
        results/haplotype_caller/contig_lists/contigs_

NUM_LISTS=$(ls results/haplotype_caller/contig_lists/ | wc -l)

CONDA="conda activate envs/gatk"


for INT in $(seq -w 1 $NUM_LISTS); do
    while read -r SAMPLE; do

        BAM_IN="results/map_reads/"$SAMPLE"_processed.bam"
        REF_IN="data/schHae2_wmito.fa"
        INT_IN="results/haplotype_caller/contig_lists/contigs_"$INT".list"
        VCF_OUT="results/haplotype_caller/"$SAMPLE"-"$INT".hc.vcf"

        CMD="bin/gatk-4.2.0.0/gatk HaplotypeCaller \
            --input $BAM_IN \
            --output $VCF_OUT \
            -reference $REF_IN \
            --emit-ref-confidence GVCF \
            -L $INT_IN"

        echo $CONDA; $CMD  | $QSUB -pe smp 2 -N scan_hc_"$SAMPLE"-"$INT" -o results/logs/hc_"$SAMPLE"-"$INT".log

    done < samples.list
done


conda activate envs/gatk

#merge ##### on high mem machine
while read -r SAMPLE; do

    ls results/haplotype_caller/"$SAMPLE"-*.hc.vcf >results/haplotype_caller/$SAMPLE.list


    bin/gatk-4.2.0.0/gatk --java-options "-Xmx128g" \
        MergeVcfs \
        --MAX_RECORDS_IN_RAM 20000000 \
        -I results/haplotype_caller/$SAMPLE.list \
        -O results/haplotype_caller/$SAMPLE.hc.vcf

done < samples.list


