conda activate scan-01-process_genome

mkdir ~/sch_hae_scan/results/vs_sb

cd ~/sch_hae_scan/results/vs_sb

BOV_REF="GCA_944470425.1"

conda activate datasets
datasets download genome accession ${BOV_REF} --filename ${BOV_REF}.zip
conda deactivate

unzip ${BOV_REF}.zip

cut -f1 -d" " ncbi_dataset/data/GCA_944470425.1/GCA_944470425.1_tdSchBovi1.1_genomic.fna >GCA_944470425.1_tdSchBovi1.1_genomic.fna

BOV_FAS="/master/nplatt/sch_hae_scan/results/vs_sb/GCA_944470425.1_tdSchBovi1.1_genomic.fna"

bwa index ${BOV_FAS}
samtools faidx ${BOV_FAS} 


QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -p -1023" 

for SAMPLE in $(); do
    CMD="conda run --live-stream -n scan-02-filter_map_genotype bin/gatk-4.2.0.0/gatk --java-options \"-Xmx36g -Xms12g\" GenomicsDBImport \
        -V results/haplotype_caller/hc.list \
        --genomicsdb-workspace-path /master/nplatt/sch_hae_scan/results/gdbimport/${CONTIG} \
        -L ${CONTIG} \
        --reader-threads 12 \
        --batch-size 12"

    #echo $CMD
    echo $CMD  | $QSUB -pe smp 12 -N gdbi_${CONTIG} -o gdbi_${CONTIG}.log

    CMD="conda run --live-stream -n scan-02-filter_map_genotype bin/gatk-4.2.0.0/gatk --java-options \"-Xmx32g -Xms18g\" GenotypeGVCFs \
                -R ${REF} \
                -V gendb:///master/nplatt/sch_hae_scan/results/gdbimport/${CONTIG} \
                -O results/genotype/${CONTIG}.vcf
                -new-qual "

    #echo $CMD
    echo $CMD  | $QSUB -pe smp 12 -N geno_${CONTIG} -o geno_${CONTIG}.log -hold_jid gdbi_${CONTIG}

done



rule pe_bbmap:
    input:
        r1_read_fq = results_dir + "/filtered_reads/{id}_filtered_R1.fq.gz",
        r2_read_fq = results_dir + "/filtered_reads/{id}_filtered_R2.fq.gz",
        ref         = genome
    output:
        pe_sam   = results_dir + "/vs_sb/{id}_pe.sam"
    threads:
        12
    log:
        logs_dir + "/pe_bbmap#{id}"
    conda:
        envs_dir + "/bbmap.yml"
    shell:
        """
        bbmap.sh \
            -ref={input.ref} \
            -nodisk \
            -in={input.r1_read_fq} \
            -in2={input.r2_read_fq} \
            vslow \
            -threads={threads} \
            ambig=toss \
            interleaved=false \
            -Xmx24g \
            -eoom \
            -out={output.pe_sam} \
            minid=0.8
        """

rule se_bbmap:
    input:
        rx_read_fq = results_dir + "/filtered_reads/{id}_filtered_RX.fq.gz",
        ref   = genome
    output:
        se_sam   = results_dir + "/vs_sb/{id}_se.sam"
    threads:
        12
    log:
        logs_dir + "/se_bbmap#{id}"
    conda:
        envs_dir + "/bbmap.yml"
    shell:
        """
        bbmap.sh \
            -ref={input.ref} \
            -nodisk \
            -in={input.rx_read_fq} \
            vslow \
            -threads={threads} \
            ambig=toss \
            -Xmx24g \
            -eoom \
            -out={output.se_sam} \
            minid=0.8
        """