REF_FAS="$SCRATCH/sch_hae_scan/data/ShV3_scaffolds.fa"
RESULTS_DIR="$SCRATCH/sch_hae_scan/results"

#mkdir $RESULTS_DIR/filter_genotypes
cd $RESULTS_DIR/filter_genotypes


conda activate vcftools

bgzip -c $RESULTS_DIR/genotype/merged_unfiltered.vcf >$RESULTS_DIR/genotype/merged_unfiltered.vcf.gz
tabix -p vcf $RESULTS_DIR/genotype/merged_unfiltered.vcf.gz


$SCRATCH/sch_hae_scan/bin/gatk-4.2.0.0/gatk VariantFiltration \
    -R $REF_FAS \
    -V $RESULTS_DIR/genotype/merged_unfiltered.vcf \
    -O $RESULTS_DIR/filter_genotypes/merged_softfiltered.vcf \
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

$SCRATCH/sch_hae_scan/bin/gatk-4.2.0.0/gatk SelectVariants \
	-V $RESULTS_DIR/filter_genotypes/merged_softfiltered.vcf \
	--exclude-filtered \
	-O $RESULTS_DIR/filter_genotypes/merged_hardfilt.vcf \
	-R $REF_FAS

conda deactivate

#get bialleleic snps with min gq and depth
conda activate envs/vcftools

vcftools \
	--vcf merged_hardfilt.vcf \
	--minGQ 15 \
	--minDP 12 \
	--min-alleles 2 \
	--max-alleles 2 \
	--remove-indels \
	--recode \
	--recode-INFO-all \
	--stdout \
	>merged_bisnps_min_gq_depth.vcf

#remove sites gtd in less than 50% of samples
vcftools \
	--vcf $RESULTS_DIR/filter_genotypes/merged_bisnps_min_gq_depth.vcf \
	--missing-site \
	 --stdout \
	>$RESULTS_DIR/filter_genotypes/missing_per_site.tbl

cat $RESULTS_DIR/filter_genotypes/missing_per_site.tbl \
	| sed 1d \
	| awk '{if ($6<=0.5) print $1"\t"$2}' \
	>$RESULTS_DIR/filter_genotypes/high_freq_gt_sites.list

vcftools \
	--vcf $RESULTS_DIR/filter_genotypes/merged_bisnps_min_gq_depth.vcf \
	--positions $RESULTS_DIR/filter_genotypes/high_freq_gt_sites.list \
	--recode \
	--recode-INFO-all \
	--stdout \
	>$RESULTS_DIR/filter_genotypes/high_freq_gt_sites.vcf

#remove individuals genotyped at less than 50% of sites
conda activate envs/vcftools

 vcftools \
	--vcf $RESULTS_DIR/filter_genotypes/high_freq_gt_sites.vcf \
	--missing-indv \
	--stdout \
	>$RESULTS_DIR/filter_genotypes/indv_gt_freq.tbl

 cat $RESULTS_DIR/filter_genotypes/indv_gt_freq.tbl \
	| awk '{if ($5<=0.5) print $1}' \
	| sed 1d \
	>$RESULTS_DIR/filter_genotypes/indvs_to_keep.list 

vcftools \
	--vcf $RESULTS_DIR/filter_genotypes/high_freq_gt_sites.vcf \
	--keep $RESULTS_DIR/filter_genotypes/indvs_to_keep.list   \
	--recode \
	--recode-INFO-all \
	--stdout \
	>$RESULTS_DIR/filter_genotypes/indv_and_site_filt.vcf

conda deactivate

#annotate remaining snps
conda activate envs/bcftools

bcftools annotate \
	--set-id +'%CHROM\:%POS' \
	$RESULTS_DIR/filter_genotypes/indv_and_site_filt.vcf \
	>$RESULTS_DIR/filter_genotypes/annotated_snps.vcf

conda deactivate

#get maf05s
conda activate envs/vcftools

vcftools \
	--vcf $RESULTS_DIR/filter_genotypes/annotated_snps.vcf \
	--maf 0.05 \
	--recode \
	--recode-INFO-all \
	--stdout \
	>$RESULTS_DIR/filter_genotypes/maf05.vcf

conda deactivate

#ld filtering
conda activate envs/plink

plink \
	--vcf $RESULTS_DIR/filter_genotypes/maf05.vcf \
	--allow-extra-chr \
	--double-id \
	--indep-pairwise 25 5 0.20 \
	--out $RESULTS_DIR/filter_genotypes/maf05_ld_filtered

conda deactivate

conda activate vcftools

vcftools \
	--vcf $RESULTS_DIR/filter_genotypes/maf05.vcf \
	--exclude $RESULTS_DIR/filter_genotypes/maf05_ld_filtered.prune.out \
	--recode \
	--recode-INFO-all \
	--stdout \
	>$RESULTS_DIR/filter_genotypes/maf05_ld_filtered.vcf

conda deactivate


#####get only ingroups (haem, sp, and bovis)
##manually edit the $RESULTS_DIR/filter_genotypes/indvs_to_keep.list 
#nano $RESULTS_DIR/filter_genotypes/indvs_to_keep.list 

conda activate envs/vcftools
vcftools \
	--vcf $RESULTS_DIR/filter_genotypes/annotated_snps.vcf \
	--keep $RESULTS_DIR/filter_genotypes/ingroup_indvs_to_keep.list   \
	--recode \
	--recode-INFO-all \
	--stdout \
	>$RESULTS_DIR/filter_genotypes/ingroup_indv_and_site_filt.vcf

#get maf05s
vcftools \
	--vcf $RESULTS_DIR/filter_genotypes/ingroup_indv_and_site_filt.vcf \
	--maf 0.05 \
	--recode \
	--recode-INFO-all \
	--stdout \
	>$RESULTS_DIR/filter_genotypes/ingroup_maf05.vcf

conda deactivate

#ld filtering
conda activate envs/plink

plink \
	--vcf $RESULTS_DIR/filter_genotypes/ingroup_maf05.vcf \
	--allow-extra-chr \
	--double-id \
	--indep-pairwise 25 5 0.20 \
	--out $RESULTS_DIR/filter_genotypes/ingroup_maf05_ld_filtered

conda deactivate

conda activate vcftools

vcftools \
	--vcf $RESULTS_DIR/filter_genotypes/ingroup_maf05.vcf \
	--exclude $RESULTS_DIR/filter_genotypes/ingroup_maf05_ld_filtered.prune.out \
	--recode \
	--recode-INFO-all \
	--stdout \
	>$RESULTS_DIR/filter_genotypes/ingroup_maf05_ld_filtered.vcf

conda deactivate

time sh -c "dd if=/dev/zero of=ddfile bs=8k count=25000 && sync"; rm ddfile