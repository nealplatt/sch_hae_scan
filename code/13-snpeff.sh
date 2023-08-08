conda activate scan_snpeff

cd /master/nplatt/sch_hae_scan/results/snpeff

java -jar snpEff.jar build -noCheckCds -noCheckProtein -gff3 -v GCF_000699445.3
java -jar -Xmx64g snpEff.jar \
    ann \
    -csvStats snpeff_stats.csv \
    -i vcf \
    -o vcf \
    -htmlStats \
    -no-downstream \
    -no-upstream \
    -fastaProt snpeff_prots.fas \
    GCF_000699445.3 \
    ../phasing/beagle.vcf \
    >snpeff_beagle.vcf


vcftools \
    --vcf snpeff_beagle.vcf \
    --chr NC_067199.1 \
    --from-bp 28392842 \
    --to-bp 29196100 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >peak.vcf
    
    
        vcftools \
        --vcf ../../filter_genotypes/mappable_maf05.vcf \
        --counts \
        --snps out_targs.list \
        --indv intercalatum_drcongo_ERR119613 \
        --indv margrebowiei_zambia_ERR310940 \
        --indv matthei_zambia_ERR103051 \
        --max-mac 0 \
        --stdout \
        >mac0_3outgroup_gts