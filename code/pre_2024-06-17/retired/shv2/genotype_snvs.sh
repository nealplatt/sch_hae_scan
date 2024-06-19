#!/usr/bin/bash

cd /master/nplatt/sch_hae_scan

mkdir results/gdbimport results/genotype

conda activate envs/bedtools
bedtools makewindows \
    -w 1000000 \
    -g data/schHae2_wmito.fa.fai \
    | awk '{print $1":"$2+1"-"$3}' \
    >results/genotype/1mb_windows.list

conda deactivate

QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -p -1023" 
CONDA="conda activate /master/nplatt/sch_hae_scan/envs/gatk"

while read -r WINDOW; do

    SAFE_WINDOW=$(echo $WINDOW | sed 's/:/_/g')
    
    #mkdir sandbox/$SAFE_WINDOW
    
    CMD="bin/gatk-4.2.0.0/gatk --java-options \"-Xmx8g -Xms4g\" GenomicsDBImport \
        -V /master/nplatt/sch_hae_scan/results/haplotype_caller/hc_vcf.list \
        --genomicsdb-workspace-path /master/nplatt/sch_hae_scan/results/gdbimport/$SAFE_WINDOW \
        -L $WINDOW \
        --reader-threads 12 \
        --batch-size 12 \
        --tmp-dir sandbox/$SAFE_WINDOW"
    
    #echo "$CONDA; $CMD"  | $QSUB -pe smp 12 -N scan_gdbimport_$SAFE_WINDOW -o results/logs/gdbimport_$SAFE_WINDOW.log
    
    CMD="bin/gatk-4.2.0.0/gatk GenotypeGVCFs \
                -R /master/nplatt/sch_hae_scan/data/schHae2_wmito.fa \
                -V gendb:///master/nplatt/sch_hae_scan/results/gdbimport/$SAFE_WINDOW \
                -new-qual \
                -O /master/nplatt/sch_hae_scan/results/genotype/$SAFE_WINDOW.vcf"
                
    echo "$CONDA; $CMD"  | $QSUB -pe smp 3 -N scan_genotype_$SAFE_WINDOW -o results/logs/genotype_$SAFE_WINDOW.log -hold_jid scan_gdbimport_$SAFE_WINDOW 


done <results/genotype/1mb_windows.list


#merge into a single unfiltered vcf

conda activate /master/nplatt/sch_hae_scan/envs/gatk

ls /master/nplatt/sch_hae_scan/results/genotype/*vcf \
    >/master/nplatt/sch_hae_scan/results/genotype/genotyped_vcfs.list
        
bin/gatk-4.2.0.0/gatk --java-options "-Xmx120g" MergeVcfs \
    --MAX_RECORDS_IN_RAM 50000000 \
    -I /master/nplatt/sch_hae_scan/results/genotype/genotyped_vcfs.list \
    -O /master/nplatt/sch_hae_scan/results/genotype/merged_unfiltered.vcf

NW_023366109.1_2000001-3000000


for DIR in NW_023366109.1_1-1000000; do

    rm -r /master/nplatt/sch_hae_scan/results/gdbimport/$DIR
    bin/gatk-4.2.0.0/gatk --java-options "-Xmx200g -Xms20g" GenomicsDBImport \
        -V /master/nplatt/sch_hae_scan/results/haplotype_caller/hc_vcf.list \
        --genomicsdb-workspace-path /master/nplatt/sch_hae_scan/results/gdbimport/$DIR \
        -L $(echo $DIR | sed 's/_/:/2') \
        --reader-threads 48 \
        --batch-size 12

    bin/gatk-4.2.0.0/gatk GenotypeGVCFs \
        -R /master/nplatt/sch_hae_scan/data/schHae2_wmito.fa \
        -V gendb:///master/nplatt/sch_hae_scan/results/gdbimport/$DIR \
        -new-qual \
        -O /master/nplatt/sch_hae_scan/results/genotype/$DIR.vcf \
        >$DIR.log 2>&1 &

done



for DIR in NW_023366109.1_1-1000000 \
        NW_023366512.1_4000001-5000000 \
        NW_023366598.1_3000001-4000000 \
        NW_023366282.1_1-4965 \
        NW_023366581.1_1-1291; do

    rm -r /master/nplatt/sch_hae_nigeria/results/gdbimport/$DIR
    bin/gatk-4.2.0.0/gatk --java-options "-Xmx200g -Xms20g" GenomicsDBImport \
        -V /master/nplatt/sch_hae_nigeria/results/haplotype_caller/hc_vcf.list \
        --genomicsdb-workspace-path /master/nplatt/sch_hae_nigeria/results/gdbimport/$DIR \
        -L $(echo $DIR | sed 's/_/:/2') \
        --reader-threads 48 \
        --batch-size 12

    bin/gatk-4.2.0.0/gatk GenotypeGVCFs \
        -R /master/nplatt/sch_hae_nigeria/data/schHae2_wmito.fa \
        -V gendb:///master/nplatt/sch_hae_nigeria/results/gdbimport/$DIR \
        -new-qual \
        -O /master/nplatt/sch_hae_nigeria/results/genotype/$DIR.vcf \
        >$DIR.log 2>&1 &

done


