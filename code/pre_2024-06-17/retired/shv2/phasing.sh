#!/bin/bash

PROJ_DIR=/master/nplatt/sch_hae_scan
RESULTS_DIR=$PROJ_DIR/results
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -p -1023" 

conda activate $PROJ_DIR/basic

cd $PROJ_DIR

mkdir $RESULTS_DIR/phasing 

grep -v "#" $RESULTS_DIR/filter_genotypes/annotated_snps.vcf \
    | awk '{printf "%s\t%s\t%.6f\t%s\n", $1, $3, $2/287000, $2}' \
    >$RESULTS_DIR/phasing/auto.map

#make contig maps
for CONTIG in $(cut -f1 $PROJ_DIR/data/schHae2_wmito.fa.fai); do

    #make genetic map
    #mkdir $RESULTS_DIR/phasing/$CONTIG
    grep $CONTIG $RESULTS_DIR/phasing/auto.map >$RESULTS_DIR/phasing/$CONTIG/$CONTIG.map

    #extract autosome specific vcf for the samples of interest
    VCF_CONDA="conda activate $PROJ_DIR/envs/vcftools"
    VCF_CMD="vcftools \
        --vcf $RESULTS_DIR/filter_genotypes/annotated_snps.vcf \
        --chr $CONTIG \
        --recode \
        --stdout \
        >$RESULTS_DIR/phasing/$CONTIG/$CONTIG.vcf"
    
    BEAGLE_CONDA="conda activate $PROJ_DIR/envs/beagle"
    BEAGLE_CMD="java -Xmx4g -jar $PROJ_DIR/envs/beagle/share/beagle-5.1_24Aug19.3e8-1/beagle.jar \
        gt=$RESULTS_DIR/filter_genotypes/annotated_snps.vcf \
        chrom=$CONTIG \
        out=$RESULTS_DIR/phasing/$CONTIG/"$CONTIG"_beagle \
        map=$RESULTS_DIR/phasing/$CONTIG/$CONTIG.map \
        nthreads=4 \
        window=60 \
        overlap=20 \
        iterations=250"
                
    #echo "$VCF_CONDA; $VCF_CMD" | $QSUB -pe smp 1 -N vcf_$CONTIG -o results/logs/vcf_$CONTIG.log
    echo "$BEAGLE_CONDA; $BEAGLE_CMD" | $QSUB -pe smp 4 -N beagle_$CONTIG -o results/logs/beagle_$CONTIG.log -hold_jid vcf_$CONTIG 
done

#wait untill all contigs are phased
grep ERROR $RESULTS_DIR/logs/beagle_N*.log \
    | cut -f1 -d":" \
    | sed 's/beagle_\(.*\).log/\1/' \
    | rev \
    | cut -f1 -d"/" \
    | rev \
    >$RESULTS_DIR/phasing/failed_phasing_runs.list

for CONTIG in $(cat $RESULTS_DIR/phasing/failed_phasing_runs.list); do

    BEAGLE_CONDA="conda activate $PROJ_DIR/envs/beagle"
    rm results/logs/beagle_$CONTIG.log
    rm -r $RESULTS_DIR/phasing/$CONTIG/"$CONTIG"_beagle*
    
    BEAGLE_CMD="java -Xmx12g -jar $PROJ_DIR/envs/beagle/share/beagle-5.1_24Aug19.3e8-1/beagle.jar \
        gt=$RESULTS_DIR/filter_genotypes/annotated_snps.vcf \
        chrom=$CONTIG \
        out=$RESULTS_DIR/phasing/$CONTIG/"$CONTIG"_beagle \
        map=$RESULTS_DIR/phasing/$CONTIG/$CONTIG.map \
        nthreads=12 \
        window=60 \
        overlap=20 \
        iterations=250"
                
    echo "$BEAGLE_CONDA; $BEAGLE_CMD" | qsub -V -cwd -S /bin/bash -q high_mem.q -j y -p -1023 -pe smp 12 -N rbeagle_$CONTIG -o results/logs/beagle_$CONTIG.log
done

gunzip $RESULTS_DIR/phasing/$CONTIG"_beagle"*vcf.gz

conda deactivate

#create one file of phased snps (for all autosomes)
conda activate $PROJ_DIR/vcftools


vcfcat \
    $(ls $RESULTS_DIR/phasing/*/*_beagle.vcf) \
    >$RESULTS_DIR/phasing/beagle.vcf

vcftools \
    --vcf $RESULTS_DIR/phasing/beagle.vcf \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >$RESULTS_DIR/phasing/beagle_maf05.vcf

conda deactivate

