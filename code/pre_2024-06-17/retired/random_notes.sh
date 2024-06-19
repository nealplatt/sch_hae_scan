cd $SCRATCH/sch_hae_scan/results/haplotype_caller

rsync --dry-run -aPh nplatt@heron.txbiomedgenetics.org:/master/nplatt/sch_hae_scan/results/haplotype_caller/ .

cd gdbimport

ls ../haplotype_caller/*.vcf >../haplotype_caller/hc.list

for contig in contigs
    bin/gatk-4.2.0.0/gatk --java-options "-Xmx48g -Xms24g" GenomicsDBImport \
        -V ../haplotype_caller/hc.list \
        --genomicsdb-workspace-path ./$CONTIG \
        -L $CONTIG \
        --reader-threads 64 \
        --batch-size 12



for VCF in $(ls *.vcf.gz); do 
    CMD="conda activate scan-02-filter_map_genotype\n\n../../bin/gatk-4.2.0.0/gatk IndexFeatureFile -I $VCF"
    echo -e $CMD >"$VCF".sh

    qsub -V -cwd -S /bin/bash -q all.q -j y -N "$VCF" -o "$VCF".log -pe smp 12 "$VCF".sh
done



cp Sb_NG_au_2.13.hc.vcf ../rerun/
cp Sb_NG_au_2.6.hc.vcf ../rerun/
cp Sb_NG_be_1.3.hc.vcf ../rerun/
cp Sb_NG_ak_2.1.hc.vcf ../rerun/
cp Sb_NG_au_2.10.hc.vcf ../rerun/
cp Sb_NG_be_1.10.hc.vcf ../rerun/
cp Sb_NG_au_2.5.hc.vcf ../rerun/
cp Sb_NG_be_1.5.hc.vcf ../rerun/
cp Sb_NG_au_1.2.hc.vcf ../rerun/

igvtools index Sb_NG_be_1.3.hc.vcf.gz


#Sb_NG_au_2.13

cd /master/nplatt/sch_hae_scan/results/haplotype_caller/rerun

for SAMPLE in Sb_NG_au_2.6 \
              Sb_NG_be_1.3 \
              Sb_NG_ak_2.1 \
              Sb_NG_au_2.10 \
              Sb_NG_be_1.10 \
              Sb_NG_au_2.5 \
              Sb_NG_be_1.5 \
              Sb_NG_au_1.2; do

    echo $SAMPLE
    #cp /master/nplatt/sch_hae_scan/results/haplotype_caller/preserve/$SAMPLE.hc.vcf .

    conda activate bcftools
    bcftools convert --threads 24 -Ob ../preserve/$SAMPLE.hc.vcf >$SAMPLE.hc.bcf
    bcftools index --threads 24 $SAMPLE.hc.bcf
    bcftools concat -a -d all -o $SAMPLE.hc.vcf.gz --threads 24 -O z $SAMPLE.hc.bcf
    tabix -fp vcf $SAMPLE.hc.vcf.gz
    conda deactivate

    conda activate igvtools
    igvtools index $SAMPLE.hc.vcf.gz
    conda deactivate

    cp $SAMPLE.hc.vcf.gz /master/nplatt/sch_hae_scan/results/haplotype_caller/
    cp $SAMPLE.hc.vcf.gz.idx /master/nplatt/sch_hae_scan/results/haplotype_caller/

    rm $SAMPLE.hc.vcf.gz.tbi 
    rm $SAMPLE.hc.bcf 
    rm $SAMPLE.hc.bcf.csi 
    rm $SAMPLE.hc.vcf

done


txRo@t@n2016her
---------------


Sb_NG_au_2.13
Sb_NG_au_2.6
Sb_NG_be_1.3
Sb_NG_ak_2.1
Sb_NG_au_2.10
Sb_NG_be_1.10
Sb_NG_au_2.5
Sb_NG_be_1.5
Sb_NG_au_1.2



gatk-launch SelectVariants \
-V gendb://test_genomicsdb -O extract.vcf.gz


for VCF in $(ls *.vcf.gz); do

    CMD="zcat $VCF \
        | cut -f1,2 \
        | grep -v \"#\" \
        | sort \
        | uniq -c \
        | awk -v vcf=$VCF '{ if (\$1>1) print VCF,\$0}' \
        >$VCF.dups"

    echo $CMD>$VCF.sh

    qsub -V -cwd -S /bin/bash -q all.q -j y -N "$VCF" -o "$VCF".log -pe smp 12 "$VCF".sh
done


done >duplicates.txt

asfd


HiC_scaffold_3 21513684
HiC_scaffold_4 7811394
HiC_scaffold_5 16768191
HiC_scaffold_5 26372642
HiC_scaffold_6 13674250
HiC_scaffold_7 19766702


 sed -i '/HiC_scaffold_1\t33771284/d' sbo_uganda_apac13_35.hc.vcf
gzip sbo_uganda_apac13_35.hc.vcf &


sed -i '/HiC_scaffold_3\t21513684/d' sha_pemba_uwandani_7.hc.vcf
echo 1
sed -i '/HiC_scaffold_4\t7811394/d' sha_pemba_uwandani_7.hc.vcf
echo 2
sed -i '/HiC_scaffold_5\t16768191/d' sha_pemba_uwandani_7.hc.vcf
echo 3
sed -i '/HiC_scaffold_5\t26372642/d' sha_pemba_uwandani_7.hc.vcf
echo 4
sed -i '/HiC_scaffold_6\t13674250/d' sha_pemba_uwandani_7.hc.vcf
echo 5
sed -i '/HiC_scaffold_7\t19766702/d' sha_pemba_uwandani_7.hc.vcf

wc -l sha_pemba_uwandani_7.hc.vcf
zcat sha_pemba_uwandani_7.hc.vcf.gz | wc -l


gzip sha_pemba_uwandani_7.hc.vcf.tmp
mv sha_pemba_uwandani_7.hc.vcf.tmp.gz sha_pemba_uwandani_7.hc.vcf.gz

tabix -fp vcf $SAMPLE.hc.vcf.gz


sbo_uganda_apac13_35.hc.vcf

for SAMPLE in sha_pemba_uwandani_7; do
    bcftools convert --threads 24 -Ob $SAMPLE.hc.vcf >$SAMPLE.hc.bcf
    bcftools index --threads 24 $SAMPLE.hc.bcf
    bcftools concat -a -d all -o $SAMPLE.hc.vcf.gz --threads 24 -O z $SAMPLE.hc.bcf
    tabix -fp vcf $SAMPLE.hc.vcf.gz
done


cd $SCRATCH/sch_hae_scan/results/gdbimport

ls $SCRATCH/sch_hae_scan/results/haplotype_caller/*.gz >hc.list

for CHR in $(seq 9 563); do
    CONTIG=HiC_scaffold_"$CHR"
    echo $CONTIG
    $SCRATCH/sch_hae_scan/gatk-4.2.0.0/gatk --java-options "-Xmx64g -Xms24g" GenomicsDBImport \
        -V hc.list \
        --genomicsdb-workspace-path $SCRATCH/sch_hae_scan/results/gdbimport/$CONTIG \
        -L $CONTIG \
        --reader-threads 64 \
        > $CONTIG.log 2>&1

done

cd $SCRATCH/sch_hae_scan/results/genotype

for CHR in $(seq 1 8); do
    CONTIG=HiC_scaffold_"$CHR"
    
    CMD=$SCRATCH"/sch_hae_scan/gatk-4.2.0.0/gatk --java-options \"-Xmx64g -Xms16g\" GenotypeGVCFs \
        -R "$SCRATCH"/sch_hae_scan/data/ShV3_scaffolds.fa \
        -V gendb://"$SCRATCH"/sch_hae_scan/results/gdbimport/"$CONTIG" \
        -O "$CONTIG".vcf.gz"

    cat header >$CONTIG.sbatch.sh
    echo $CMD >>$CONTIG.sbatch.sh

    sbatch $CONTIG.sbatch.sh

done


for CHR in $(seq 9 563); do
    CONTIG=HiC_scaffold_"$CHR"

    echo $CONTIG
    
    $SCRATCH/sch_hae_scan/gatk-4.2.0.0/gatk --java-options "-Xmx64g -Xms16g" GenotypeGVCFs \
        -R $SCRATCH/sch_hae_scan/data/ShV3_scaffolds.fa \
        -V gendb://$SCRATCH/sch_hae_scan/results/gdbimport/$CONTIG \
        -O $CONTIG.vcf.gz \
        > $CONTIG.log 2>&1

done


#!/bin/bash
#----------------------------------------------------
#SBATCH -o myjob.o%j                # Name of stdout output file
#SBATCH -e myjob.e%j                # Name of stderr error file
#SBATCH -p normal                   # Queue (partition) name
#SBATCH -N 1                        # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                        # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 48:00:00                 # Run time (hh:mm:ss)
#SBATCH --mail-user=rplatt@txbiomed.org
#SBATCH --mail-type=all             # Send email at begin and end of job
#SBATCH -A ascaris_introgressio     # Allocation name (req'd if you have more than 1)


CONTIG=HiC_scaffold_8
$SCRATCH/sch_hae_scan/gatk-4.2.0.0/gatk --java-options "-Xmx80g -Xms24g" GenomicsDBImport \
    -V hc.list \
    --genomicsdb-workspace-path $SCRATCH/sch_hae_scan/results/gdbimport/nobatch/$CONTIG/ \
    -L $CONTIG \
    --reader-threads 64


CONTIG=HiC_scaffold_8
$SCRATCH/sch_hae_scan/gatk-4.2.0.0/gatk --java-options "-Xmx64g -Xms24g" GenomicsDBImport \
    -V $SCRATCH/sch_hae_scan/results/gdbimport/hc.list \
    --genomicsdb-workspace-path $SCRATCH/sch_hae_scan/results/gdbimport/nobatch/$CONTIG \
    -L $CONTIG \
    --reader-threads 64

CONTIG=HiC_scaffold_8
$SCRATCH/sch_hae_scan/gatk-4.2.0.0/gatk --java-options "-Xmx64g -Xms24g" GenomicsDBImport \
    -V ../hc.list \
    --genomicsdb-workspace-path ./$CONTIG \
    -L $CONTIG \
    --reader-threads 64


############################################################################################
cd $SCRATCH/sch_hae_scan/results/gdbimport

../../gatk-4.2.0.0/gatk SplitIntervals \
    --output intervals/ \
    -R ../../data/ShV3_scaffolds.fa \
    --subdivision-mode INTERVAL_SUBDIVISION \
    --scatter-count 100

ls $SCRATCH/sch_hae_scan/results/haplotype_caller/*.gz >hc.list

#make header
#!/bin/bash
#----------------------------------------------------
#SBATCH -o myjob.o%j                # Name of stdout output file
#SBATCH -e myjob.e%j                # Name of stderr error file
#SBATCH -p normal                   # Queue (partition) name
#SBATCH -N 1                        # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                        # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 48:00:00                 # Run time (hh:mm:ss)
#SBATCH --mail-user=rplatt@txbiomed.org
#SBATCH --mail-type=all             # Send email at begin and end of job
#SBATCH -A ascaris_introgressio     # Allocation name (req'd if you have more than 1)


for INT in $(seq -w 0049 0098); do
    
    INT_FILE=/scratch/08286/nplatt/sch_hae_scan/results/gdbimport_intervals/"$INT"-scattered.interval_list 

    CMD=$SCRATCH"/sch_hae_scan/gatk-4.2.0.0/gatk --java-options \"-Xmx64g -Xms16g\" GenomicsDBImport \
        -V hc.list \
        --genomicsdb-workspace-path 0$INT \
        -L $INT_FILE \
        --reader-threads 4 \
        --batch-size 50 \
        >"$INT".log 2>&1"

    cat header >$INT.sbatch.sh
    echo $CMD >>$INT.sbatch.sh

    sbatch $INT.sbatch.sh

done


for INT in $(seq -w 0 0048); do
    sbatch $INT.sbatch.sh
done

for INT in $(seq -w 0049 0098); do
    sbatch $INT.sbatch.sh
done

cd $SCRATCH/sch_hae_scan/results/genotype

for INT in $(seq -w 0050 0098); do
    #INT_FILE=/scratch/08286/nplatt/sch_hae_scan/results/gdbimport_intervals/"$INT"-scattered.interval_list 

    
    #CMD=$SCRATCH"/sch_hae_scan/gatk-4.2.0.0/gatk --java-options \"-Xmx64g -Xms16g\" GenotypeGVCFs \
    #    -R "$SCRATCH"/sch_hae_scan/data/ShV3_scaffolds.fa \
    #    -V gendb://"$SCRATCH"/sch_hae_scan/results/gdbimport/0"$INT" \
    #    -O "$INT".vcf.gz \
    #    >"$INT".log 2>&1"

    #cat header >$INT.sbatch.sh
    #echo $CMD >>$INT.sbatch.sh

    sbatch $INT.sbatch.sh
done

#running
for FQ in $(ls *_1.fastq | head -n 35); do
    gzip -c $FQ >$FQ.gz &
done

#done
for FQ in $(ls *_1.fastq | tail -n 33); do
    gzip -c $FQ >$FQ.gz &
done

#running
for FQ in $(ls *_2.fastq | head -n 35); do
    gzip -c $FQ >$FQ.gz &
done

#running
for FQ in $(ls *_2.fastq | tail -n 33); do
    gzip -c $FQ >$FQ.gz &
done