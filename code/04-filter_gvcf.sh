cd /master/nplatt/sch_hae_scan

mkdir results/filter_genotypes
cd results/filter_genotypes/

REF_FAS="/master/nplatt/sch_hae_scan/data/GCF_000699445.3_UoM_Shae.V3_genomic.fna"
RESULTS_DIR="/master/nplatt/sch_hae_scan/results"

conda activate scan-02-filter_map_genotype

~/sch_hae_scan/bin/gatk-4.2.0.0/gatk VariantFiltration \
    -R ${REF_FAS} \
    -V ${RESULTS_DIR}/genotype/merge_unfiltered.vcf.gz \
    -O ./merged_softfiltered.vcf.gz \
    --filter-expression "QD < 2.0" \
    --filter-name "qd_lt_2" \
    --filter-expression "MQ < 30.0" \
    --filter-name "mq_gt_30" \
    --filter-expression "FS > 60.0" \
    --filter-name "fs_lt_60" \
    --filter-expression "SOR > 3.0" \
    --filter-name "SOR_lt_3" \
    --filter-expression "MQRankSum < -12.5" \
    --filter-name "MQRankSum_lt_-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" \
    --filter-name "ReadPosRankSum_lt_-8"

#hard filter snps
~/sch_hae_scan/bin/gatk-4.2.0.0/gatk SelectVariants \
    -V merged_softfiltered.vcf.gz \
    --exclude-filtered \
    --select-type-to-include SNP \
    -O merged_hardfilt.vcf.gz \
    -R ${REF_FAS}


#get bialleleic snps with min gq and depth
vcftools \
    --gzvcf merged_hardfilt.vcf.gz \
    --minGQ 15 \
    --minDP 10 \
    --min-alleles 2 \
    --max-alleles 2 \
    --remove-indels \
    --recode \
    --recode-INFO-all \
    --stdout \
    >merged_bisnps_min_gq_depth.vcf.gz

#remove sites gtd in less than 50% of samples
vcftools \
    --vcf merged_bisnps_min_gq_depth.vcf \
    --missing-site \
    --stdout \
    >missing_per_site.tbl

cat missing_per_site.tbl \
    | sed 1d \
    | awk '{if ($6<=0.5) print $0}' \
    >high_freq_gt_sites.list

vcftools \
    --vcf merged_bisnps_min_gq_depth.vcf \
    --positions high_freq_gt_sites.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >high_freq_gt_sites.vcf

#remove individuals genotyped at less than 50% of sites
#conda activate envs/vcftools

 vcftools \
    --vcf high_freq_gt_sites.vcf \
    --missing-indv \
    --stdout \
    >indv_gt_freq.tbl

cat indv_gt_freq.tbl \
    | awk '{if ($5<=0.5) print $1}' \
    | sed 1d \
    >indvs_to_keep.list 

vcftools \
    --vcf high_freq_gt_sites.vcf \
    --keep indvs_to_keep.list   \
    --recode \
    --recode-INFO-all \
    --stdout \
    >indv_and_site_filt.vcf


#annotate
bcftools annotate \
    --set-id +'%CHROM\:%POS' \
    indv_and_site_filt.vcf \
    >annotated_snps.vcf

../../bin/gatk-4.2.0.0/gatk SortVcf --java-options "-Xmx256g -Xms24g" \
      I=annotated_snps.vcf \
      O=sorted_annotated_snps.vcf

# run the find_mappable_regionss.sh script to get a white list set of regions.
sed 's/ Schistosoma haematobium chromosome .*sequence//' ../mappable_regions/mappable_regions.bed >shv3_mappable_regions.bed

# vcftools \
#     --vcf sorted_annotated_snps.vcf \
#     --bed shv3_mappable_regions.bed \
#     --recode \
#     --recode-INFO-all \
#     --stdout \
#     >mappable_snvs.vcf
    
# #on to phasing