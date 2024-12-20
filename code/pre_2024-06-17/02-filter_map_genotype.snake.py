import os
import pandas as pd
from Bio import SeqIO
from snakemake.utils import min_version

#using
# - snakemake v6.9.1
# - using conda v22.9.0

#ex:
#snakemake \
# --printshellcmds \
# --cluster 'qsub -V -cwd -j y -S /bin/bash -pe smp {threads} -q all.q -o {log} ' \
# --jobs 200 \
# --latency-wait 200 \
# --keep-going \
# --rerun-incomplete \
# --snake code/02-filter_map_genotype.snake \
# --use-conda \
# --jobname snk.{name}.{wildcards.id}.jid{jobid}

##### set minimum snakemake version #####
min_version("6.9.1")

#set main project dir and work from there
proj_dir    = "/master/nplatt/sch_hae_scan"
data_dir    = "{}/data".format(proj_dir)
results_dir = "{}/results".format(proj_dir)
envs_dir    = "{}/envs".format(proj_dir)
logs_dir    = "{}/logs".format(results_dir)

#set ref genome
genome = "{}/GCF_000699445.3_UoM_Shae.V3_genomic.fna".format(data_dir)

#get sample info
samples_df  = pd.read_csv("{}/samples.csv".format(data_dir))
samples     = list(samples_df["sample_id"])

#get contigs
contigs = []
with open("{}.fai".format(genome), 'r') as in_fai:
    for entry in in_fai:
        contig = entry.rstrip().split("\t")[0]
        contigs.append(contig)

localrules: 
    all, 

rule all:
    input:
        expand("{dir}/seq_data/{id}_{read}.fq.gz", dir  = data_dir, 
                                                   id   = samples, 
                                                   read = ["R1", "R2"] ),
        expand("{dir}/filtered_reads/{id}_filtered_{read}.fq.gz", dir  = results_dir,
                                                                  id   = samples,
                                                                  read = ["R1", "R2"]),
        # expand("{dir}/mapped_reads/{id}_pe.sam", dir  = results_dir,
        #                                          id   = samples),
        # expand("{dir}/mapped_reads/{id}_se.sam", dir  = results_dir,
        #                                          id   = samples),
        # expand("{dir}/mapped_reads/{id}_processed{ext}", dir = results_dir, 
        #                                               id  = samples, 
        #                                               ext = [".bam", ".bam.bai"] ),
        # expand("{dir}/haplotype_caller/{sample}.hc.vcf",  dir    = results_dir,
        #                                                   sample = samples),
        # expand("{dir}/gdbimport/{contig}/{files}", dir    = results_dir, 
        #                                            contig = contigs,
        #                                            files  = [ "callset.json", 
        #                                                       "__tiledb_workspace.tdb", 
        #                                                       "vcfheader.vcf", 
        #                                                       "vidmap.json" ]),
        # expand("{dir}/genotype/{contig}.vcf", dir     = results_dir,
        #                                       contig  = contigs),
        # results_dir + "/genotype/merged_unfiltered.vcf",
        # expand("{dir}/mapped_reads/{id}_pe_flagstat.txt", dir  = results_dir,
        #                                                   id   = samples),

rule filter_reads:
    input:
        r1           = data_dir + "/seq_data/{id}_R1.fq.gz",
        r2           = data_dir + "/seq_data/{id}_R2.fq.gz",
        adapter_file = data_dir + "/adapters.fas"
    output:
        r1_pe = results_dir + "/filtered_reads/{id}_filtered_R1.fq.gz",
        r2_pe = results_dir + "/filtered_reads/{id}_filtered_R2.fq.gz",
        rx_se = results_dir + "/filtered_reads/{id}_filtered_RX.fq.gz",
        r1_se = results_dir + "/filtered_reads/{id}_filtered_SE_R1.fq.gz",
        r2_se = results_dir + "/filtered_reads/{id}_filtered_SE_R2.fq.gz"
    threads:
        12
    log:
        logs_dir + "/filter_reads#{id}.log"
    conda:
        envs_dir + "/trimmomatic.yml"
    shell:
        """
        trimmomatic \
            PE \
            -threads {threads} \
            -phred33 \
            {input.r1} \
            {input.r2} \
            {output.r1_pe} \
            {output.r1_se} \
            {output.r2_pe} \
            {output.r2_se} \
            LEADING:10 \
            TRAILING:10 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 \
            ILLUMINACLIP:{input.adapter_file}:2:30:10:1:true

        zcat {output.r1_se} {output.r2_se} | gzip >{output.rx_se} 
        """

rule pe_bbmap:
    input:
        r1_read_fq = results_dir + "/filtered_reads/{id}_filtered_R1.fq.gz",
        r2_read_fq = results_dir + "/filtered_reads/{id}_filtered_R2.fq.gz",
        ref         = genome
    output:
        pe_sam   = results_dir + "/mapped_reads/{id}_pe.sam"
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
        se_sam   = results_dir + "/mapped_reads/{id}_se.sam"
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

rule sort_merge_bbmap:
    input:
        pe_sam   = results_dir + "/mapped_reads/{id}_pe.sam",
        se_sam   = results_dir + "/mapped_reads/{id}_se.sam"
    output:
        pe_bam_sorted = temp(results_dir + "/mapped_reads/{id}_sorted_pe.bam"),
        se_bam_sorted = temp(results_dir + "/mapped_reads/{id}_sorted_se.bam"),
        merged_bam    = temp(results_dir + "/mapped_reads/{id}_merged.bam")
    threads:
        12
    log:
        logs_dir + "/sort_merge_bbmap#{id}"
    conda:
        envs_dir + "/samtools.yml"
    shell:
      """
        samtools view -Sb {input.se_sam} | samtools sort -o {output.se_bam_sorted}
        samtools view -Sb {input.pe_sam} | samtools sort -o {output.pe_bam_sorted}
            
        samtools merge \
            {output.merged_bam} \
            {output.pe_bam_sorted} \
            {output.se_bam_sorted} \
        """

# rule add_readgroups:
#     input:
#         merged_bam = results_dir + "/mapped_reads/{id}_merged.bam"
#     output:
#         rg_bam = temp(results_dir + "/mapped_reads/{id}_rg.bam")
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
#         rg_sortred_bam = temp(results_dir + "/mapped_reads/{id}_rg_sorted.bam")
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
#         metrics = temp(results_dir + "/mapped_reads/{id}_dupmetrics.log"),
#         bam     = results_dir + "/mapped_reads/{id}_processed.bam",
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
#         bai = results_dir + "/mapped_reads/{id}_processed.bam.bai"
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

rule flagstat:
    input:
        pe_sam   = results_dir + "/mapped_reads/{id}_pe.sam"
    output:
        txt = results_dir + "/mapped_reads/{id}_pe_flagstat.txt"
    threads:
        1
    log:
        logs_dir + "/flagstat#{id}"
    conda:
        envs_dir + "/samtools.yml"
    shell:
        """
        samtools flagstat {input.pe_sam} >{output.txt}
        """

# rule haplotype_caller:
#     input:
#         bai = results_dir + "/mapped_reads/{id}_processed.bam.bai",
#         bam = results_dir + "/mapped_reads/{id}_processed.bam",
#         ref = genome
#     output:
#         vcf = results_dir + "/haplotype_caller/{id}.hc.vcf"
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
#         vcfs = expand("{dir}/haplotype_caller/{id}.hc.vcf", dir    = results_dir, 
#                                                             id  = samples) 
#     output:
#         vcf_list = results_dir + "/haplotype_caller/hc.list",
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
# rule mk_gdbimport_dir:
#     output:
#         results_dir + "/gdbimport"
#     threads:
#         1
#     log:
#         logs_dir + "/mk_gdbimport_dir"
#     shell:
#         """
#         mkdir {output}
#         """

# rule gdbimport:
#     input:
#         rules.mk_gdbimport_dir.output,
#         vcf_list = results_dir + "/haplotype_caller/hc.list",
#     output:
#         expand("{dir}/gdbimport/{{contig}}/callset.json", dir = results_dir), 
#                                                      #files =  [ "callset.json", 
#                                                      #           "__tiledb_workspace.tdb", 
#                                                      #           "vcfheader.vcf", 
#                                                      #           "vidmap.json" ]),
#     params:
#         gdb_dir = results_dir + "/gdbimport/{contig}"
#     threads:
#         12
#     log:
#         logs_dir + "/gdbimport_{contig}"
#     conda:
#         envs_dir + "/scan-02-filter_map_genotype.yml"
#     shell:
#         """
#         bin/gatk-4.2.0.0/gatk --java-options \"-Xmx36g -Xms12g\" GenomicsDBImport \
#                 -V {input.vcf_list} \
#                 --genomicsdb-workspace-path {params.gdb_dir} \
#                 -L {wildcards.contig} \
#                 --reader-threads {threads} \
#                 --batch-size 12
#         """

# rule genotype:
#     input:
#         results_dir + "/gdbimport/{contig}/callset.json",
#         vcf_list = results_dir + "/haplotype_caller/hc.list",
#         ref = genome
#     output:
#         vcf = results_dir + "/genotype/{contig}.vcf"
#     params:
#         gdb_dir = results_dir + "/gdbimport/{contig}"
#     threads:
#         12
#     log:
#         logs_dir + "/genotype_{contig}"
#     conda:
#         envs_dir + "/scan-02-filter_map_genotype.yml"
#     shell:
#         """
#         bin/gatk-4.2.0.0/gatk --java-options \"-Xmx32g -Xms18g\" GenotypeGVCFs \
#                 -R {input.ref} \
#                 -V gendb://{params.gdb_dir} \
#                 -O {output.vcf}
#                 -new-qual \
#         """    

# rule merge_gvcf:
#     input:
#         vcfs = expand(results_dir + "/genotype/{contig}.vcf", contig=contigs),
#         ref  = genome
#     output:
#         vcf_list = results_dir + "/genotype/hc.list",
#         vcf = results_dir + "/genotype/merged_unfiltered.vcf"
#     threads:
#         12
#     log:
#         logs_dir + "/merge_gvcf"
#     conda:
#         envs_dir + "/scan-02-filter_map_genotype.yml"
#     shell:
#         """
#         ls {input.vcfs} >{output.vcf_list}

#         bin/gatk-4.2.0.0/gatk --java-options \"-Xmx120g\" MergeVcfs \
#             --MAX_RECORDS_IN_RAM 50000000 \
#             -I {output.vcf_list} \
#             -O {output.vcf}
#         """   
