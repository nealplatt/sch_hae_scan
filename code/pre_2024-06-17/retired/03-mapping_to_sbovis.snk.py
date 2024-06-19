import os
import pandas as pd
from Bio import SeqIO
from snakemake.utils import min_version

#using
# - snakemake v6.9.1
# - using conda v22.9.0

# ex:
# snakemake \
# --printshellcmds \
# --cluster 'qsub -V -cwd -j y -S /bin/bash -pe smp {threads} -q all.q -o {log} ' \
# --jobs 200 \
# --latency-wait 200 \
# --keep-going \
# --rerun-incomplete \
# --snake code/mapping_to_sbovis.snk.py \
# --use-conda \
# --jobname snk.{name}.{wildcards.contig}.jid{jobid}

##### set minimum snakemake version #####
min_version("6.9.1")

#set main project dir and work from there
proj_dir    = "/master/nplatt/sch_hae_scan"
data_dir    = "{}/data".format(proj_dir)
results_dir = "{}/results".format(proj_dir)
envs_dir    = "{}/envs".format(proj_dir)
logs_dir    = "{}/logs".format(results_dir)

#set ref genome
genome = "{}/vs_sb/GCA_944470425.1_tdSchBovi1.1_genomic.fna".format(results_dir)

#get sample info
samples_df  = pd.read_csv("{}/samples.csv".format(data_dir))
samples     = list(samples_df["sample_id"])

#get contigs
contigs = []
with open("{}.fai".format(genome), 'r') as in_fai:
    for entry in in_fai:
        contig = entry.rstrip().split("\t")[0]
        contigs.append(contig)


# conda activate scan-01-process_genome

# mkdir ~/sch_hae_scan/results/vs_sb

# cd ~/sch_hae_scan/results/vs_sb

# BOV_REF="GCA_944470425.1"

# conda activate datasets
# datasets download genome accession ${BOV_REF} --filename ${BOV_REF}.zip
# conda deactivate

# unzip ${BOV_REF}.zip

# cut -f1 -d" " ncbi_dataset/data/GCA_944470425.1/GCA_944470425.1_tdSchBovi1.1_genomic.fna >GCA_944470425.1_tdSchBovi1.1_genomic.fna

# BOV_FAS="/master/nplatt/sch_hae_scan/results/vs_sb/GCA_944470425.1_tdSchBovi1.1_genomic.fna"

# bwa index ${BOV_FAS}
# samtools faidx ${BOV_FAS} 


localrules: 
    all, 

rule all:
    input:
        # expand("{dir}/filtered_reads/{id}_filtered_{read}.fq.gz", dir  = results_dir,
        #                                                           id   = samples,
        #                                                           read = ["R1", "R2"]),
        # expand("{dir}/vs_sb/{id}_pe.sam", dir  = results_dir,
        #                                   id   = samples),
        # expand("{dir}/vs_sb/{id}_se.sam", dir  = results_dir,
        #                                   id   = samples),
        # expand("{dir}/vs_sb/{id}_processed{ext}", dir = results_dir, 
        #                                           id  = samples, 
        #                                           ext = [".bam", ".bam.bai"] ),
        # expand("{dir}/vs_sb/{id}_pe_flagstat.txt", dir  = results_dir,
        #                                            id   = samples),
        # expand("{dir}/vs_sb/hc/{id}.hc.vcf", dir  = results_dir,
        #                                   id   = samples),
        expand("{dir}/vs_sb/gdbimport/{contig}/{files}", dir    = results_dir, 
                                               contig = contigs,
                                               files  = [ "callset.json", 
                                                          "__tiledb_workspace.tdb", 
                                                          "vcfheader.vcf", 
                                                          "vidmap.json" ]),
        expand("{dir}/vs_sb/genotype/{contig}.vcf", dir     = results_dir,
                                           contig  = contigs),
        results_dir + "/vs_sb/merged_unfiltered.vcf",
        # expand("{dir}/mapped_reads/{id}_pe_flagstat.txt", dir  = results_dir,
        #                                                   id   = samples),

# rule pe_bbmap:
#     input:
#         r1_read_fq = results_dir + "/filtered_reads/{id}_filtered_R1.fq.gz",
#         r2_read_fq = results_dir + "/filtered_reads/{id}_filtered_R2.fq.gz",
#         ref         = genome
#     output:
#         pe_sam   = results_dir + "/vs_sb/{id}_pe.sam"
#     threads:
#         12
#     log:
#         logs_dir + "/pe_bbmap#{id}"
#     conda:
#         envs_dir + "/bbmap.yml"
#     shell:
#         """
#         bbmap.sh \
#             -ref={input.ref} \
#             -nodisk \
#             -in={input.r1_read_fq} \
#             -in2={input.r2_read_fq} \
#             vslow \
#             -threads={threads} \
#             ambig=toss \
#             interleaved=false \
#             -Xmx24g \
#             -eoom \
#             -out={output.pe_sam} \
#             minid=0.8
#         """

# rule se_bbmap:
#     input:
#         rx_read_fq = results_dir + "/filtered_reads/{id}_filtered_RX.fq.gz",
#         ref   = genome
#     output:
#         se_sam   = results_dir + "/vs_sb/{id}_se.sam"
#     threads:
#         12
#     log:
#         logs_dir + "/se_bbmap#{id}"
#     conda:
#         envs_dir + "/bbmap.yml"
#     shell:
#         """
#         bbmap.sh \
#             -ref={input.ref} \
#             -nodisk \
#             -in={input.rx_read_fq} \
#             vslow \
#             -threads={threads} \
#             ambig=toss \
#             -Xmx24g \
#             -eoom \
#             -out={output.se_sam} \
#             minid=0.8
#         """

# rule sort_merge_bbmap:
#     input:
#         pe_sam   = results_dir + "/vs_sb/{id}_pe.sam",
#         se_sam   = results_dir + "/vs_sb/{id}_se.sam"
#     output:
#         pe_bam_sorted = temp(results_dir + "/vs_sb/{id}_sorted_pe.bam"),
#         se_bam_sorted = temp(results_dir + "/vs_sb/{id}_sorted_se.bam"),
#         merged_bam    = temp(results_dir + "/vs_sb/{id}_merged.bam")
#     threads:
#         12
#     log:
#         logs_dir + "/sort_merge_bbmap#{id}"
#     conda:
#         envs_dir + "/samtools.yml"
#     shell:
#       """
#         samtools view -Sb {input.se_sam} | samtools sort -o {output.se_bam_sorted}
#         samtools view -Sb {input.pe_sam} | samtools sort -o {output.pe_bam_sorted}
            
#         samtools merge \
#             {output.merged_bam} \
#             {output.pe_bam_sorted} \
#             {output.se_bam_sorted} \
#         """

# rule add_readgroups:
#     input:
#         merged_bam = results_dir + "/vs_sb/{id}_merged.bam"
#     output:
#         rg_bam = temp(results_dir + "/vs_sb/{id}_rg.bam")
#     threads:
#         1
#     log:
#         logs_dir + "/add_readgroups#{id}"
#     conda:
#         envs_dir + "/scan-02-filter_map_genotype.yml"
#     shell:
#         """
#         bin/gatk-4.2.0.0/gatk AddOrReplaceReadGroups \
#             --INPUT {input.merged_bam} \
#             --OUTPUT {output.rg_bam} \
#             --RGPU unk \
#             --RGLB library1 \
#             --RGPL illumina \
#             --RGSM {wildcards.id} \
#             --RGID {wildcards.id}
#         """

# rule sort_bam_post_readgroup:
#     input:
#         rg_bam = rules.add_readgroups.output.rg_bam,
#     output:
#         rg_sortred_bam = temp(results_dir + "/vs_sb/{id}_rg_sorted.bam")
#     threads:
#         12
#     log:
#         logs_dir + "/sort_bam_post_readgroup#{id}"
#     conda:
#         envs_dir + "/samtools.yml"
#     shell:
#         """
#         samtools sort -o {output.rg_sortred_bam} {input.rg_bam}
#         """

# # #cheater code here to create a file for coverage estimates

# rule mark_duplicates_in_bam:
#     input:
#         rg_sortred_bam = rules.sort_bam_post_readgroup.output.rg_sortred_bam
#     output:
#         metrics = temp(results_dir + "/vs_sb/{id}_dupmetrics.log"),
#         bam     = results_dir + "/vs_sb/{id}_processed.bam",
#     threads:
#         2
#     log:
#         logs_dir + "/mark_duplicates_in_bam#{id}"
#     conda:
#         envs_dir + "/scan-02-filter_map_genotype.yml"
#     shell:
#         """
#         bin/gatk-4.2.0.0/gatk MarkDuplicates \
#             --INPUT {input.rg_sortred_bam} \
#             --OUTPUT {output.bam} \
#             --METRICS_FILE {output.metrics} \
#             --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900 \
#             --ASSUME_SORT_ORDER coordinate
#         """

# rule index_bam:
#     input:
#         bam = rules.mark_duplicates_in_bam.output.bam
#     output:
#         bai = results_dir + "/vs_sb/{id}_processed.bam.bai"
#     threads:
#         1
#     log:
#         logs_dir + "/index_bam#{id}"
#     conda:
#         envs_dir + "/samtools.yml"
#     shell:
#         """
#         samtools index {input.bam}
#         """

# rule flagstat:
#     input:
#         pe_sam   = results_dir + "/vs_sb/{id}_pe.sam"
#     output:
#         txt = results_dir + "/vs_sb/{id}_pe_flagstat.txt"
#     threads:
#         1
#     log:
#         logs_dir + "/flagstat#{id}"
#     conda:
#         envs_dir + "/samtools.yml"
#     shell:
#         """
#         samtools flagstat {input.pe_sam} >{output.txt}
#         """

# rule haplotype_caller:
#     input:
#         bai = results_dir + "/vs_sb/{id}_processed.bam.bai",
#         bam = results_dir + "/vs_sb/{id}_processed.bam",
#         ref = genome
#     output:
#         vcf = results_dir + "/vs_sb/{id}.hc.vcf"
#     threads:
#         12
#     log:
#         logs_dir + "/haplotype_caller_{id}"
#     conda:
#         envs_dir + "/scan-02-filter_map_genotype.yml"
#     shell:
#         """
#         bin/gatk-4.2.0.0/gatk --java-options \"-Xmx36g\" HaplotypeCaller \
#             --input {input.bam} \
#             --output {output.vcf} \
#             -reference {input.ref} \
#             --emit-ref-confidence GVCF \
#             --native-pair-hmm-threads {threads}
#         """

# rule list_hc_vcfs:
#     input:
#         vcfs = expand("{dir}/vs_sb/{id}.hc.vcf", dir = results_dir, 
#                                                  id  = samples) 
#     output:
#         vcf_list = results_dir + "/vs_sb/hc.list",
#     threads:
#         12
#     log:
#         logs_dir + "/list_hc_vcfs"
#     conda:
#         envs_dir + "/scan-02-filter_map_genotype.yml"
#     shell:
#         """
#         ls {input.vcfs}>{output.vcf_list}
#         """


rule gdbimport:
    input:
        vcf_list = results_dir + "/vs_sb/hc.list",
    output:
        expand("{dir}/vs_sb/gdbimport/{{contig}}/{files}", dir = results_dir, 
                                                     files =  [ "callset.json", 
                                                                "__tiledb_workspace.tdb", 
                                                                "vcfheader.vcf", 
                                                                "vidmap.json" ]),
        gdbi_dir = directory(results_dir + "/vs_sb/gdbimport/{contig}")
    threads:
        12
    log:
        logs_dir + "/gdbimport_{contig}.log"
    conda:
        envs_dir + "/scan-02-filter_map_genotype.yml"
    shell:
        """
        bin/gatk-4.2.0.0/gatk --java-options \"-Xmx36g -Xms12g\" GenomicsDBImport \
                -V {input.vcf_list} \
                --genomicsdb-workspace-path {output.gdbi_dir} \
                -L {wildcards.contig} \
                --reader-threads {threads} \
                --batch-size 12
        """

rule genotype:
    input:
        gdbi_dir = results_dir + "/vs_sb/gdbimport/{contig}",
        vcf_list = results_dir + "/vs_sb/hc.list",
        ref = genome
    output:
        vcf = results_dir + "/vs_sb/genotype/{contig}.vcf"
    threads:
        12
    log:
        logs_dir + "/vs_sb_{contig}"
    conda:
        envs_dir + "/scan-02-filter_map_genotype.yml"
    shell:
        """
        bin/gatk-4.2.0.0/gatk --java-options \"-Xmx32g -Xms18g\" GenotypeGVCFs \
                -R {input.ref} \
                -V gendb://{input.gdbi_dir} \
                -O {output.vcf}
                -new-qual \
        """    

rule merge_gvcf:
    input:
        vcfs = expand(results_dir + "/vs_sb/genotype/{contig}.vcf", contig=contigs),
        ref  = genome
    output:
        vcf_list = results_dir + "/vs_sb/genotype/hc.list",
        vcf = results_dir + "/vs_sb/merged_unfiltered.vcf"
    threads:
        12
    log:
        logs_dir + "/merge_gvcf"
    conda:
        envs_dir + "/scan-02-filter_map_genotype.yml"
    shell:
        """
        ls {input.vcfs} >{output.vcf_list}

        bin/gatk-4.2.0.0/gatk --java-options \"-Xmx120g\" MergeVcfs \
            --MAX_RECORDS_IN_RAM 50000000 \
            -I {output.vcf_list} \
            -O {output.vcf}
        """   


# samtools view -Sb -h -F 4 f1sbsh_unk_SRR7743803_pe.sam   >f1_mapped.bam
# samtools sort -n -o f1_mapped_sorted.bam f1_mapped.bam
# bedtools bamtofastq -i f1_mapped_sorted.bam -fq f1_mapped_to_sb_R1.fq -fq2 f1_mapped_to_sb_R2.fq

# samtools view -Sb -h -f 4 f1sbsh_unk_SRR7743803_pe.sam >f1_unmapped.bam
# samtools sort -n -o f1_unmapped_sorted.bam f1_unmapped.bam
# bedtools bamtofastq -i f1_unmapped_sorted.bam -fq f1_unmapped_to_sb_R1.fq -fq2 f1_unmapped_to_sb_R2.fq


# bbmap.sh \
#     -ref=/master/nplatt/sch_hae_scan/data/GCF_000699445.3_UoM_Shae.V3_genomic.fna \
#     -in=f1_mapped_to_sb_R1.fq \
#     -in2=f1_mapped_to_sb_R2.fq \
#     vslow \
#     -threads=24 \
#     ambig=toss \
#     interleaved=false \
#     -Xmx24g \
#     -eoom \
#     -out=f1_mapped_to_sb_mapped_back_to_sh.sam \
#     minid=0.8


# bbmap.sh \
#     -ref=/master/nplatt/sch_hae_scan/data/GCF_000699445.3_UoM_Shae.V3_genomic.fna \
#     -in=f1_unmapped_to_sb_R1.fq \
#     -in2=f1_unmapped_to_sb_R2.fq \
#     vslow \
#     -threads=24 \
#     ambig=toss \
#     interleaved=false \
#     -Xmx24g \
#     -eoom \
#     -out=f1_unmapped_to_sb_mapped_back_to_sh.sam \
#     minid=0.8

# samtools flagstat f1_mapped_to_sb_mapped_back_to_sh.sam >f1_mapped_to_sb_mapped_back_to_sh.flagstat
# samtools flagstat f1_unmapped_to_sb_mapped_back_to_sh.sam >f1_unmapped_to_sb_mapped_back_to_sh.flagstat

# #get mapping rates
# cd /master/nplatt/sch_hae_scan/results/vs_sb

# for TXT in $(ls *_pe_flagstat.txt); do

#     SAMPLE=$(echo ${TXT} | sed 's/_pe_flagstat.txt//')
#     SH=$(grep -m1 mapped ../mapped_reads/${TXT} | cut -f2 -d"(" | cut -f1 -d" ")
#     SB=$(grep -m1 mapped ${TXT}| cut -f2 -d"(" | cut -f1 -d" ")
#     echo $SAMPLE,$SH,$SB
# done



# for BAM in $(ls *_processed.bam); do
    
#     QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 24 -N ${BAM} -o ${BAM}tocram.log"
#     CRAM=$(echo $BAM | sed 's/.bam/.cram/')
#     REF="/master/nplatt/sch_hae_scan/results/vs_sb/GCA_944470425.1_tdSchBovi1.1_genomic.fna"
    
#     CMD="conda run -n scan-02-filter_map_genotype samtools view --fai-reference  ${REF} --cram --output ${CRAM} ${BAM}"

#     echo $CMD | $QSUB

# done
