#convert bam to cram

REF="/master/nplatt/sch_hae_scan/results/vs_sb/GCA_944470425.1_tdSchBovi1.1_genomic.fna"
REF="master/nplatt/sch_hae_scan/data/GCF_000699445.3_UoM_Shae.V3_genomic.fna"

for BAM in $(ls *_processed.bam); do

    CRAM=$(echo $BAM | sed 's/.bam/.cram/')
    echo $BAM
    conda run -n scan-02-filter_map_genotype samtools view --fai-reference  ${REF} --cram --output ${CRAM} ${BAM}

done
