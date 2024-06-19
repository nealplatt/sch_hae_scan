conda activate scan-04-phasing

mkdir /master/nplatt/sch_hae_scan/results/phasing
cd /master/nplatt/sch_hae_scan/results/phasing

#calculate a recomb map.  1cm = 287Kb based on data from Criscione et al. (2009) avg
grep -v "#" ../filter_genotypes/sorted_annotated_snps.vcf  \
    | awk '{printf "%s\t%s\t%.6f\t%s\n", $1, $3, $2/287000, $2}' \
    >gen.map

#get contigs with ge 2 variants
cut -f1 gen.map | sort | uniq -c >chr.counts

awk '{if ($1>=100) print $2}' chr.counts >ge100_chr.counts

#make dirs for each contig/chrom
for CHR in $(cat chr.counts  | awk '{print $2}'); do
    mkdir $CHR
done

#print chrom specific recomb maps
awk '{print $0 > $1"/"$1".gen.map" }' gen.map


# cat header.sge.sh
# #$ -cwd
# #$ -V
# #$ -S /bin/bash
# #$ -N CHR
# #$ -o CHR.stdeo
# #$ -j y
# #$ -q all.q
# #$ -pe smp 48

# conda activate scan-04-phasing


#impute on each seperatley (create seperate qsub scripts)
for CHR in $(cat ge100_chr.counts); do

    cd $CHR

    sed "s/CHR/$CHR/" ../header.sge.sh >$CHR.sge.sh

    #extract autosome specific vcf for the samples of interest
    CMD="vcftools \
        --vcf /master/nplatt/sch_hae_scan/results/filter_genotypes/sorted_annotated_snps.vcf \
        --chr $CHR \
        --recode \
        --recode-INFO-all \
        --stdout \
        >$CHR.vcf"
    echo -e $CMD"\n" >>$CHR.sge.sh

    # CMD="sed '1,606s/,assembly=ShV3_scaffolds.fa//gi' $CHR.vcf >$CHR.mod.vcf"
    # echo -e $CMD"\n" >>$CHR.sge.sh

    CMD="java -Xmx64g -jar /master/nplatt/anaconda3/envs/scan-04-phasing/share/beagle-5.2_21Apr21.304-0/beagle.jar \
        gt=$CHR.vcf \
        out="$CHR"_beagle \
        map=$CHR.gen.map \
        nthreads=48 \
        window=20 \
        overlap=10 \
        iterations=1000 \
        burnin=1000 \
        >$CHR.log 2>&1"
    echo -e $CMD"\n" >>$CHR.sge.sh

    CMD="gunzip "$CHR"_beagle.vcf.gz"
    echo -e $CMD"\n" >>$CHR.sge.sh

    qsub $CHR.sge.sh
    cd ..

done

#wait for everything to phase

#create one file of phased snps (for all autosomes)
vcfcombine $(ls HiC_scaffold_*/*_beagle.vcf) >beagle.vcf


vcftools \
    --vcf beagle.vcf \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >beagle_maf05.vcf
    #41,680