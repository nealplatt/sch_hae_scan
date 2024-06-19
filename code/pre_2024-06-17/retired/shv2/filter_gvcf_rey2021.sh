cd /master/nplatt/sch_hae_scan

mkdir results/filter_genotypes_for_f1

cd results/filter_genotypes_for_f1

REF_FAS="/master/nplatt/sch_hae_scan/data/schHae2_wmito.fa"
RESULTS_DIR="/master/nplatt/sch_hae_scan/results"

conda activate /master/nplatt/sch_hae_scan/envs/gatk

/master/nplatt/sch_hae_scan/bin/gatk-4.2.0.0/gatk VariantFiltration \
    --QUIET True \
    --verbosity ERROR \
    -R $REF_FAS \
    -V $RESULTS_DIR/genotype/merged_unfiltered.vcf \
    -O merged_softfiltered_rey2021.vcf \
    --filter-expression "QD < 2.0" \
    --filter-name "qd_lt_2" \
    --filter-expression "MQ < 40.0" \
    --filter-name "mq_gt_30" \
    --filter-expression "FS > 45.0" \
    --filter-name "fs_lt_45" \
    --filter-expression "MQRankSum < -12.5" \
    --filter-name "MQRankSum_lt_-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" \
    --filter-name "ReadPosRankSum_lt_-8"


#hard filter snps
/master/nplatt/sch_hae_scan/bin/gatk-4.2.0.0/gatk SelectVariants \
    -V merged_softfiltered_rey2021.vcf \
    --exclude-filtered \
    -O merged_hardfilt_rey2021.vcf \
    -R $REF_FAS

#filter based on depth etc. -> rey dataset
conda activate /master/nplatt/sch_hae_scan/envs/vcftools

#get the good sites from the rey f1
 vcftools \
    --vcf merged_hardfilt_rey2021.vcf \
    --indv f1sbsh_unk_SRR7743803 \
    --minDP 8 \
    --remove-indels \
    --recode \
    --recode-INFO-all \
    --stdout \
    --remove-filtered-all \
    >f1_rey2021.vcf 

vcftools \
	--vcf f1_rey2021.vcf \
	--missing-site \
	 --stdout \
	>f1_rey2021_per_site.tbl

cat f1_rey2021_per_site.tbl \
	| sed 1d \
	| awk '{if ($6==0) print $1"\t"$2}' \
	>f1_rey2021_per_site_keep.list

vcftools \
    --vcf merged_hardfilt_rey2021.vcf \
    --positions f1_rey2021_per_site_keep.list \
    --minDP 8 \
    --mac 3 \
    --min-alleles 2 \
    --max-alleles 2 \
    --remove-indels \
    --recode \
    --recode-INFO-all \
    --stdout \
    --remove-filtered-all \
	>rey2021_site_filtered.vcf

vcftools \
    --vcf rey2021_site_filtered.vcf \
    --max-missing 0.5 \
    --recode \
    --recode-INFO-all \
    --stdout \
    --remove-filtered-all \
    >rey2021_site_filtered_max50p.vcf

#remove individuals genotyped at less than 50% of sites
 vcftools \
    --vcf rey2021_site_filtered_max50p.vcf \
    --missing-indv \
    --stdout \
    >indv_gt_freq.tbl

 cat indv_gt_freq.tbl \
    | awk '{if ($5<=0.5) print $1}' \
    | sed 1d \
    >indvs_to_keep.list 

#manually removed outgroups

vcftools \
    --vcf rey2021_site_filtered_max50p.vcf \
    --keep indvs_to_keep.list   \
    --recode \
    --recode-INFO-all \
    --stdout \
    >rey2021_indv_and_site_filt.vcf


conda activate  /master/nplatt/sch_hae_scan/envs/bcftools

bcftools annotate \
    --set-id +'%CHROM\:%POS' \
    rey2021_indv_and_site_filt.vcf \
    >rey2021_filtered.vcf

conda deactivate


conda activate  /master/nplatt/sch_hae_scan/envs/vcftools

vcftools \
    --vcf rey2021_filtered.vcf \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >rey2021_maf05.vcf

#ld filtering
conda activate /master/nplatt/sch_hae_scan/envs/plink

plink \
    --vcf rey2021_maf05.vcf \
    --allow-extra-chr \
    --double-id \
    --indep-pairwise 25 5 0.20 \
    --out rey2021_maf05_ld_filtered

conda deactivate

conda activate /master/nplatt/sch_hae_scan/envs/vcftools

vcftools \
    --vcf rey2021_maf05.vcf \
    --exclude rey2021_maf05_ld_filtered.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >rey2021_maf05_ld_filtered.vcf

conda deactivate


#now run PCA and admixture