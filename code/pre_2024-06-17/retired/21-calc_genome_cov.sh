QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 4 " 
REF="/master/nplatt/sch_hae_scan/data/GCF_000699445.3_UoM_Shae.V3_genomic.fna"

mkdir /master/nplatt/sch_hae_scan/results/mosdepth
cd /master/nplatt/sch_hae_scan/results/mosdepth

for CRAM in $(ls /master/nplatt/sch_hae_scan/results/mapped_reads/*.cram); do

    SAMPLE=$(basename $CRAM .cram)
                
    CMD="conda run -n mosdepth mosdepth --threads 4 --no-per-base --fasta $REF --use-median $SAMPLE $CRAM"
        
    echo $CMD | $QSUB -N md.$SAMPLE -o $SAMPLE.mosdepth.log
done


##################################################################
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 1 " 

mkdir /master/nplatt/sch_hae_scan/results/read_counts_from_fq/
cd /master/nplatt/sch_hae_scan/results/read_counts_from_fq/

for FQGZ in $(ls /master/nplatt/sch_hae_scan/data/seq_data/*R1.fq.gz); do

    SAMPLE=$(basename $FQGZ _R1.fq.gz)
    NUM=$(( $(zcat $FQGZ | wc -l) / 4 ))
                
    echo $SAMPLE,$NUM 
done>read_counts_from_fq.csv

