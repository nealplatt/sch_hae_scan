QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -p -1023" 
REF="/master/nplatt/sch_hae_scan/data/GCF_000699445.3_UoM_Shae.V3_genomic.fna"

for CONTIG in $(cut -f1 data/GCF_000699445.3_UoM_Shae.V3_genomic.fna.fai); do

    CMD="conda run --live-stream -n scan-02-filter_map_genotype \
            bin/gatk-4.2.0.0/gatk --java-options \"-Xmx32g -Xms18g\" GenotypeGVCFs \
                -R ${REF} \
                -V gendb:///master/nplatt/sch_hae_scan/results/gdbimport/${CONTIG} \
                -O results/genotype/${CONTIG}.vcf
                -new-qual "

    #echo $CMD
    echo $CMD  | $QSUB -pe smp 12 -N geno_${CONTIG} -o geno_${CONTIG}.log -hold_jid gdbi_${CONTIG}

done


#merge vcfs

#compress haplotype caller
for SAMPLE in $(ls *hc.vcf | sed 's/.hc.vcf//'); do
    CMD="conda run --live-stream -n scan-02-filter_map_genotype /master/nplatt/sch_hae_scan/bin/gatk-4.2.0.0/gatk --java-options \"-Xmx32g -Xms18g\" VcfFormatConverter INPUT=${SAMPLE}.hc.vcf OUTPUT=${SAMPLE}.hc.vcf.gz REQUIRE_INDEX=true"

    echo $CMD  | $QSUB -pe smp 4 -N bcf_${SAMPLE} -o bcf_${SAMPLE}.log

done


#merge vcfs




for VCF in $(ls *hc.vcf); do
    echo $VCF
    bgzip -c ${VCF} >${VCF}.gz
    tabix -p vcf ${VCF}.gz
done


#compress haplotype caller
for VCF in $(ls *hc.vcf); do
    CMD="conda run --live-stream -n scan-02-filter_map_genotype bgzip -c ${VCF} >${VCF}.gz"

    echo $CMD  | $QSUB -pe smp 2 -N vcfgz_${VCF} -o vcfgz_${VCF}.log

done

for VCFGZ in $(ls *hc.vcf.gz); do
    CMD="conda run --live-stream -n scan-02-filter_map_genotype tabix -p vcf ${VCFGZ}"

    echo $CMD  | $QSUB -pe smp 2 -N tabix_${VCFGZ} -o tabix_${VCFGZ}.log

    CMD="conda run --live-stream -n scan-02-filter_map_genotype /master/nplatt/sch_hae_scan/bin/gatk-4.2.0.0/gatk --java-options \"-Xmx10g -Xms2g\" IndexFeatureFile -I ${VCFGZ}"

    echo $CMD  | $QSUB -pe smp 4 -N idx_${VCFGZ} -o idx_${VCFGZ}.log

done
