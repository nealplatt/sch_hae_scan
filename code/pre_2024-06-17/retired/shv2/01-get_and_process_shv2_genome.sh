#!/usr/bin/bash

#get in proper dir and activate conda env
cd ~/sch_hae_scan
conda activate envs/get_and_process_sh_genome

GENOME_ACCESSION=GCF_000699445.2
#move into dir for sra data
cd data

ncbi-genome-download \
    --section refseq \
    --formats fasta,gff \
    --assembly-accessions $GENOME_ACCESSION \
    invertebrate

gunzip refseq/invertebrate/GCF_000699445.2/GCF_000699445.2_SchHae_2.0_genomic.fna.gz

#mask mitochondrial genome?
        #------------------------------------------------------
        # VERY IMPORTANT NOTE
        # From previous work I know that the mitochondria has been missasembled into a larger contig
        # I am masking that region here and adding the mitochondria from NCBI.  This can be 
        # verifed by doing a simple blast search and finding seeing that there is no single
        # contig that has a large Sh mitochondrial hit....but is within the expected
        # size range of a mitochondria (8-30kb). This info is stored in 
        # data/missassembled_mito.bed
        #------------------------------------------------------

echo -e "AMPZ02000201.1\t1\t20391" >missassembled_mito.bed
bedtools maskfasta \
    -fi refseq/invertebrate/GCF_000699445.2/GCF_000699445.2_SchHae_2.0_genomic.fna \
    -fo tmp_schHae2.fa \
    -bed missassembled_mito.bed

#add new mito
esearch -db nucleotide -query NC_008074.1 |  efetch -format fasta >NC_008074.1_mito.fas
cat tmp_schHae2.fa NC_008074.1_mito.fas >schHae2_wmito.fa

#remove bad mito
rm tmp_schHae2.fa

#-------------------
# Build indecies/dicts

#bwa index
bwa index schHae2_wmito.fa

#gatk seq dict
conda activte ../envs/gatk
../bin/gatk-4.2.0.0/gatk CreateSequenceDictionary -R schHae2_wmito.fa
conda deactivate

#faidx
samtools faidx schHae2_wmito.fa


