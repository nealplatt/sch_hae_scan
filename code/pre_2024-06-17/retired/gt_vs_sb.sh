
conda activate scan-02-filter_map_genotype

cd /master/nplatt/sch_hae_scan/results/vs_sb

QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 12" 

REF="/master/nplatt/sch_hae_scan/results/vs_sb/GCA_944470425.1_tdSchBovi1.1_genomic.fna"

ls $(pwd)/hc/*hc.vcf >hc.list
cut -f1 ${REF}.fai >contigs.list

mkdir gdbimport genotype

for CONTIG in $(cat contigs.list); do

    CMD="conda run --cwd . --no-capture-output -n scan-02-filter_map_genotype \
        /master/nplatt/sch_hae_scan/bin/gatk-4.2.0.0/gatk --java-options \"-Xmx36g -Xms12g\" GenomicsDBImport \
            -V /master/nplatt/sch_hae_scan/results/vs_sb/hc.list \
            --genomicsdb-workspace-path /master/nplatt/sch_hae_scan/results/vs_sb/gdbimport/${CONTIG} \
            -L ${CONTIG} \
            --reader-threads 12 \
            --batch-size 12"

    #echo $CMD  | $QSUB -N gdbi_${CONTIG} -o gdbi_${CONTIG}.log

    CMD="conda run --cwd . --no-capture-output -n scan-02-filter_map_genotype \
        /master/nplatt/sch_hae_scan/bin/gatk-4.2.0.0/gatk --java-options \"-Xmx32g -Xms18g\" GenotypeGVCFs \
                -R ${REF} \
                -V gendb:///gdbimport/${CONTIG} \
                -O /master/nplatt/sch_hae_scan/results/vs_sb/genotype/${CONTIG}.vcf
                -new-qual"

    echo $CMD  | $QSUB -N gvcf_${CONTIG} -o gvcf_${CONTIG}.log -hold_jid gdbi_${CONTIG}

done


ls /master/nplatt/sch_hae_scan/results/vs_sb/genotype/*.vcf >gvcfs.list

/master/nplatt/sch_hae_scan/bin/gatk-4.2.0.0/gatk --java-options "-Xmx120g" MergeVcfs \
    --MAX_RECORDS_IN_RAM 50000000 \
    -I /master/nplatt/sch_hae_scan/results/vs_sb/gvcfs.list \
    -O /master/nplatt/sch_hae_scan/results/vs_sb/genotype/merged_unfiltered_vs_sb.vcf

####################################################################################
mkdir filter_genotypes
cd filter_genotypes/

REF="/master/nplatt/sch_hae_scan/results/vs_sb/GCA_944470425.1_tdSchBovi1.1_genomic.fna"


conda activate scan-02-filter_map_genotype

~/sch_hae_scan/bin/gatk-4.2.0.0/gatk VariantFiltration \
    -R ${REF_FAS} \
    -V /master/nplatt/sch_hae_scan/results/vs_sb/genotype/merged_unfiltered_vs_sb.vcf.gz \
    -O ./merged_softfiltered_vs_sb.vcf.gz \
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
    -V merged_softfiltered_vs_sb.vcf.gz \
    --exclude-filtered \
    --select-type-to-include SNP \
    -O merged_hardfilt_vs_sb.vcf.gz \
    -R ${REF_FAS}


#get bialleleic snps with min gq and depth
vcftools \
    --gzvcf merged_hardfilt_vs_sb.vcf.gz \
    --minGQ 15 \
    --minDP 10 \
    --min-alleles 2 \
    --max-alleles 2 \
    --remove-indels \
    --recode \
    --recode-INFO-all \
    --stdout \
    >merged_bisnps_min_gq_depth_vs_sb.vcf.gz


#remove sites gtd in less than 50% of samples
vcftools \
    --vcf merged_bisnps_min_gq_depth_vs_sb.vcf \
    --missing-site \
    --stdout \
    >missing_per_site_vs_sb.tbl

cat missing_per_site.tbl \
    | sed 1d \
    | awk '{if ($6<=0.5) print $0}' \
    >high_freq_gt_sites_vs_sb.list

vcftools \
    --vcf merged_bisnps_min_gq_depth_vs_sb.vcf \
    --positions high_freq_gt_sites_vs_sb.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >high_freq_gt_sites_vs_sb.vcf

#remove individuals genotyped at less than 50% of sites
#conda activate envs/vcftools

 vcftools \
    --vcf high_freq_gt_sites_vs_sb.vcf \
    --missing-indv \
    --stdout \
    >indv_gt_freq_vs_sb.tbl

cat indv_gt_freq_vs_sb.tbl \
    | awk '{if ($5<=0.5) print $1}' \
    | sed 1d \
    >indvs_to_keep_vs_sb.list 

vcftools \
    --vcf high_freq_gt_sites_vs_sb.vcf \
    --keep indvs_to_keep_vs_sb.list   \
    --recode \
    --recode-INFO-all \
    --stdout \
    >indv_and_site_filt_vs_sb.vcf


#annotate
bcftools annotate \
    --set-id +'%CHROM\:%POS' \
    indv_and_site_filt_vs_sb.vcf \
    >annotated_snps_vs_sb.vcf

~/sch_hae_scan/bin/gatk-4.2.0.0/gatk SortVcf --java-options "-Xmx256g -Xms24g" \
      I=annotated_snps_vs_sb.vcf \
      O=sorted_annotated_snps_vs_sb.vcf

#get maf05s
vcftools \
    --vcf sorted_annotated_snps_vs_sb.vcf \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >maf05_vs_sb.vcf

#ld filtering
plink \
    --vcf maf05_vs_sb.vcf \
    --allow-extra-chr \
    --double-id \
    --indep-pairwise 25 5 0.20 \
    --out maf05_ld_filtered_vs_sb

#Pruning complete.  6801863 of 7231736 variants removed.
vcftools \
    --vcf maf05_vs_sb.vcf \
    --exclude maf05_ld_filtered_vs_sb.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >maf05_ld_filtered_vs_sb.vcf

#compress large vcf files
for i in merged_bisnps_min_gq_depth_vs_sb.vcf \
         high_freq_gt_sites_vs_sb.vcf \
         indv_and_site_filt_vs_sb.vcf \
         sorted_annotated_snps_vs_sb.vcf; do
    bgzip -c $i >$i.gz &
done

wait

for i in merged_bisnps_min_gq_depth_vs_sb.vcf \
         high_freq_gt_sites_vs_sb.vcf \
         indv_and_site_filt_vs_sb.vcf \
         sorted_annotated_snps_vs_sb.vcf; do
    tabix -p vcf $i.gz &
done


###############################################################################################################################
# run the find_mappable_regionss.sh script to get a white list set of regions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!sed 's/ Schistosoma haematobium chromosome .*sequence//' ../mappable_regions/mappable_regions.bed >shv3_mappable_regions.bed

vcftools \
    --vcf sorted_annotated_snps_vs_sb.vcf \
    --bed shv3_mappable_regions_vs_sb.bed \
    --recode \
    --recode-INFO-all \
    --stdout \
    >mappable_snvs_vs_sb.vcf


#get maf05s
vcftools \
    --vcf mappable_snvs_vs_sb.vcf \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >mappable_maf05_vs_sb.vcf

#ld filtering
plink \
    --vcf mappable_maf05_vs_sb.vcf \
    --allow-extra-chr \
    --double-id \
    --indep-pairwise 25 5 0.20 \
    --out mappable_maf05_ld_filtered_vs_sb

#Pruning complete.  6801863 of 7231736 variants removed.
vcftools \
    --vcf mappable_maf05_vs_sb.vcf \
    --exclude mappable_maf05_ld_filtered_vs_sb.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >mappable_maf05_ld_filtered_vs_sb.vcf


(bcftools) (11:10 AM) [nplatt@zeus filter_genotypes]$ bcftools query -l maf05_ld_filtered.vcf >samples.list
(bcftools) (11:10 AM) [nplatt@zeus filter_genotypes]$ nano samples.list
(bcftools) (11:11 AM) [nplatt@zeus filter_genotypes]$ bcftools query -l maf05_ld_filtered.vcf >samples.list^C
(bcftools) (11:11 AM) [nplatt@zeus filter_genotypes]$ mv samples.list ingroup.list
(bcftools) (11:11 AM) [nplatt@zeus filter_genotypes]$ bcftools query -l maf05_ld_filtered.vcf >samples.list


vcftools \
    --vcf mappable_maf05_vs_sb.vcf \
    --keep ingroup.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >mappable_maf05_ld_filtered_vs_sb.vcf