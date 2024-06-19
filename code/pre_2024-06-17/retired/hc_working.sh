# curassoni_senegal_ERR119623 - re-mapping
# bovis_keyna_SRR13579878     - running hc
# bovis_ethiopia_SRR13579874  - running hc

PROJ_DIR="/master/nplatt/sch_hae_scan"
SAMPLE=curassoni_senegal_ERR119623

for i in $(seq -w 0 141); do
        R1_FQ=reads/curassoni_senegal_ERR119623_1m-"$i"_filtered_R1.fq
        R2_FQ=reads/curassoni_senegal_ERR119623_1m-"$i"_filtered_R2.fq
        REF=/master/nplatt/sch_hae_scan/data/ShV3_scaffolds.fa
        SAM=curassoni_senegal_ERR119623_1m-"$i".sam

        CMD="conda activate scan-02-filter_map_genotype; \n\n bbmap.sh \
            -ref=$REF \
            -in="$R1_FQ" \
            -in2="$R2_FQ" \
            vslow \
            -threads=12 \
            ambig=toss \
            interleaved=false \
            -Xmx8g \
            -eoom \
            -out="$SAM" \
            minid=0.8"

        echo -e $CMD >bbmap-"$i".sh

        qsub -V -cwd -S /bin/bash -j y -N bbmap-"$i" -o bbmap-"$i".log -pe smp 12 bbmap-"$i".sh

done


#remove unmapped and sort 1m
while [ $(qstat | grep bbmap | wc -l) -gt 0 ]; do
    echo -ne "."
    sleep 1m
done

 for i in $(seq -w 0 141); do
        R1_FQ=reads/curassoni_senegal_ERR119623_1m-"$i"_filtered_R1.fq
        R2_FQ=reads/curassoni_senegal_ERR119623_1m-"$i"_filtered_R2.fq
        REF=/master/nplatt/sch_hae_scan/data/ShV3_scaffolds.fa
        SAM=curassoni_senegal_ERR119623_1m-"$i".sam
        BAM=curassoni_senegal_ERR119623_1m-"$i".bam

        samtools view -F 4 -Sb $SAM | samtools sort -o $BAM

done

#samtools view -Sb curassoni_senegal_ERR119623_RX.sam | samtools sort -o curassoni_senegal_ERR119623_RX.bam

#merge
samtools merge -f \
    curassoni_senegal_ERR119623_merged.bam \
    $(ls curassoni_senegal_ERR119623_1m-???.bam) \
    curassoni_senegal_ERR119623_RX.bam

#replace read groups
/master/nplatt/sch_hae_scan/bin/gatk-4.2.0.0/gatk AddOrReplaceReadGroups \
    --INPUT curassoni_senegal_ERR119623_merged.bam \
    --OUTPUT curassoni_senegal_ERR119623_merged_RG.bam \
    --RGPU unk \
    --RGLB library1 \
    --RGPL illumina \
    --RGSM curassoni_senegal_ERR119623 \
    --RGID curassoni_senegal_ERR119623

#re-sort
samtools sort -o curassoni_senegal_ERR119623_merged_RG_sorted.bam curassoni_senegal_ERR119623_merged_RG.bam

#mark dups
/master/nplatt/sch_hae_scan/bin/gatk-4.2.0.0/gatk MarkDuplicates \
    --INPUT curassoni_senegal_ERR119623_merged_RG_sorted.bam \
    --OUTPUT curassoni_senegal_ERR119623_processed.bam \
    --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900 \
    --ASSUME_SORT_ORDER coordinate \
    --METRICS_FILE curassoni_senegal_ERR119623.metrics

cp curassoni_senegal_ERR119623_processed.bam /master/nplatt/sch_hae_scan/results/mapped_reads/

#haplotype caller
SAMPLE=curassoni_senegal_ERR119623
for i in $(seq -w 0000 0098); do

    PROJ_DIR="/master/nplatt/sch_hae_scan"

    CMD="conda activate scan-02-filter_map_genotype\n\n"$PROJ_DIR"/bin/gatk-4.2.0.0/gatk --java-options \"-Xmx8g\" HaplotypeCaller \
        --input "$PROJ_DIR"/results/mapped_reads/"$SAMPLE"_processed.bam \
        --output "$PROJ_DIR"/scratch/"$i"-"$SAMPLE".vcf \
        -reference "$PROJ_DIR"/data/ShV3_scaffolds.fa \
        --emit-ref-confidence GVCF \
        --native-pair-hmm-threads 12 \
        -L "$PROJ_DIR"/scratch/"$i"-scattered.interval_list"

    echo -e $CMD >gatk-"$i"-"$SAMPLE".sh

    qsub -V -cwd -S /bin/bash -q all.q -j y -N gatk-"$i"-"$SAMPLE" -o gatk-"$i"-"$SAMPLE".log -pe smp 12 gatk-"$i"-"$SAMPLE".sh
done


#get list of vcfs
##########################################################################################
SAMPLE=curassoni_senegal_ERR119623


ls "$PROJ_DIR"/scratch/00*-"$SAMPLE".vcf >"$PROJ_DIR"/scratch/"$SAMPLE"_vcfs.list

"$PROJ_DIR"/bin/gatk-4.2.0.0/gatk MergeVcfs \
    -I "$PROJ_DIR"/scratch/"$SAMPLE"_vcfs.list \
    -O "$PROJ_DIR"/scratch/"$SAMPLE"-unsorted.vcf

"$PROJ_DIR"/bin/gatk-4.2.0.0/gatk SortVcf \
      I="$PROJ_DIR"/scratch/"$SAMPLE"-unsorted.vcf \
      O="$PROJ_DIR"/scratch/"$SAMPLE".hc.vcf

mv "$PROJ_DIR"/scratch/"$SAMPLE".hc.vcf* "$PROJ_DIR"/results/haplotype_caller/
##########################################################################################


../../bin/gatk-4.2.0.0/gatk --java-options "-Xmx8g" HaplotypeCaller \
        --input curassoni_senegal_ERR119623_merged_RG_sorted.bam \
        --output test.vcf \
        -reference ../../data/ShV3_scaffolds.fa \
        --emit-ref-confidence GVCF \
        --native-pair-hmm-threads 12 \
        -L ../0000-scattered.interval_list

../../bin/gatk-4.2.0.0/gatk --java-options "-Xmx8g" HaplotypeCaller \
        --input curassoni_senegal_ERR119623_processed.bam \
        --output test.vcf \
        -reference ../../data/ShV3_scaffolds.fa \
        --emit-ref-confidence GVCF \
        --native-pair-hmm-threads 12 \
        -L ../0000-scattered.interval_list

SAMPLE=curassoni_senegal_ERR119623


#for i in $(seq -w 0000 0098); do
#for i in "0077" "0086" "0094" "0095" "0097"; do
for i in "0086" "0094" "0095"; do
for i in "0097"; do

    PROJ_DIR="/master/nplatt/sch_hae_scan"

    rm gatk-"$i"-"$SAMPLE".log
    CMD="conda activate scan-02-filter_map_genotype\n\n"$PROJ_DIR"/bin/gatk-4.2.0.0/gatk --java-options \"-Xmx8g\" HaplotypeCaller \
        --input "$PROJ_DIR"/scratch/curassoni/curassoni_senegal_ERR119623_merged_RG_sorted.bam \
        --output "$PROJ_DIR"/scratch/"$i"-"$SAMPLE".vcf \
        -reference "$PROJ_DIR"/data/ShV3_scaffolds.fa \
        --emit-ref-confidence GVCF \
        --native-pair-hmm-threads 12 \
        -L "$PROJ_DIR"/scratch/"$i"-scattered.interval_list"

    echo -e $CMD >gatk-"$i"-"$SAMPLE".sh

    #qsub -V -cwd -S /bin/bash -q all.q -j y -N gatk-"$i"-"$SAMPLE" -o gatk-"$i"-"$SAMPLE".log -pe smp 12 gatk-"$i"-"$SAMPLE".sh
done
###################################################



