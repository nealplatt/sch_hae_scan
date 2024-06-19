cd /master/nplatt/sch_hae_scan

mkdir results/filter_genotypes

REF_FAS="$SCRATCH/sch_hae_scan/data/ShV3_scaffolds.fa"
RESULTS_DIR="$SCRATCH/sch_hae_scan/results"

conda activate envs/gatk

$SCRATCH/sch_hae_scan/bin/gatk-4.2.0.0/gatk VariantFiltration \
    -R $REF_FAS \
    -V $RESULTS_DIR/genotype/merged_unfiltered.vcf.gz \
    -O $RESULTS_DIR/filter_genotypes/merged_softfiltered.vcf.gz \
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

bin/gatk-4.2.0.0/gatk SelectVariants \
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

conda deactivate

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
	--gzvcf high_freq_gt_sites.vcf.gz \
	--missing-indv \
	--stdout \
	>indv_gt_freq.tbl

cat indv_gt_freq.tbl \
	| awk '{if ($5<=0.5) print $1}' \
	| sed 1d \
	>indvs_to_keep.list 

vcftools \
	--gzvcf high_freq_gt_sites.vcf.gz \
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

../../bin/gatk-4.2.0.0/gatk SortVcf --java-options "-Xmx48g -Xms12g" \
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


for i in annotated_snps.vcf \
         indv_and_site_filt.vcf \
         merged_bisnps_min_gq_depth.vcf \
         sorted_annotated_snps.vcf; do
    #bgzip -c $i >$i.gz &
    tabix -p vcf $i.gz &
done

bwa aln 



################################ PHASING

cd /master/nplatt/sch_hae_scan/results/phasing
conda activate scan-03-filter_genotypes

for CONTIG in $(cut -f1 ../../data/ShV3_scaffolds.fa.fai); do
    echo $CONTIG
    mkdir $CONTIG

    cd $CONTIG
    CMD="scan-03-filter_genotypes; vcftools \
        --gzvcf ../../filter_genotypes/annotated_snps.vcf \
        --chr $CONTIG \
        --recode \
        --recode-INFO-all \
        --stdout \
        >$CONTIG.vcf"

    echo $CMD | qsub -V -cwd -S /bin/bash -j y -N vcf-"$CONTIG" -o vcf-"$CONTIG".log -pe smp 1 

    cd ..

done