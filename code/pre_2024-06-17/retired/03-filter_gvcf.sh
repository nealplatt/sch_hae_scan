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

#get maf05s
vcftools \
    --vcf sorted_annotated_snps.vcf \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >maf05.vcf

#ld filtering
plink \
    --vcf maf05.vcf \
    --allow-extra-chr \
    --double-id \
    --indep-pairwise 25 5 0.20 \
    --out maf05_ld_filtered

#Pruning complete.  6801863 of 7231736 variants removed.
vcftools \
	--vcf maf05.vcf \
	--exclude maf05_ld_filtered.prune.out \
	--recode \
	--recode-INFO-all \
	--stdout \
	>maf05_ld_filtered.vcf

#####get only ingroups (haem, sp, and bovis)
# create a text file with outgroup samples to be removed (`outgroups.list`)
bcftools query -l maf05_ld_filtered.vcf >samples.list
#nano samples.list >outgroups.list

vcftools \
    --vcf maf05_ld_filtered.vcf \
    --remove outgroup.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >sh_sb_maf05_ld_filtered.vcf


#compress large vcf files
for i in merged_bisnps_min_gq_depth.vcf \
         high_freq_gt_sites.vcf \
         indv_and_site_filt.vcf \
         sorted_annotated_snps.vcf; do
    bgzip -c $i >$i.gz &
    #tabix -p vcf $i.gz &
done


for BAM in $(ls *_processed.bam); do
    
    QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 24 -N ${BAM} -o ${BAM}tocram.log"
    CRAM=$(echo $BAM | sed 's/.bam/.cram/')
    REF="/master/nplatt/sch_hae_scan/data/GCF_000699445.3_UoM_Shae.V3_genomic.fna"
    
    CMD="conda run -n scan-02-filter_map_genotype samtools view --fai-reference  ${REF} --cram --output ${CRAM} ${BAM}"

    echo $CMD | $QSUB

done

###############################################################################################################################
# run the find_mappable_regionss.sh script to get a white list set of regions.
sed 's/ Schistosoma haematobium chromosome .*sequence//' ../mappable_regions/mappable_regions.bed >shv3_mappable_regions.bed

vcftools \
    --vcf sorted_annotated_snps.vcf \
    --bed shv3_mappable_regions.bed \
    --recode \
    --recode-INFO-all \
    --stdout \
    >mappable_snvs.vcf


#get maf05s
vcftools \
    --vcf mappable_snvs.vcf \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >mappable_maf05.vcf

#ld filtering
plink \
    --vcf mappable_maf05.vcf \
    --allow-extra-chr \
    --double-id \
    --indep-pairwise 25 5 0.20 \
    --out mappable_maf05_ld_filtered

#Pruning complete.  6801863 of 7231736 variants removed.
vcftools \
    --vcf mappable_maf05.vcf \
    --exclude mappable_maf05_ld_filtered.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >mappable_maf05_ld_filtered.vcf

#####get only ingroups (haem, sp, and bovis)
# create a text file with outgroup samples to be removed (`outgroups.list`)
bcftools query -l mappable_maf05_ld_filtered.vcf >mappable_samples.list
#nano mappable_samples.list >mappable_outgroups.list


vcftools \
    --vcf mappable_maf05_ld_filtered.vcf \
    --remove mappable_outgroups.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >mappable_sh_sb_maf05_ld_filtered.vcf
