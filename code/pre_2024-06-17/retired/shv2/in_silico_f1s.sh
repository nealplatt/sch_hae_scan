#in silico F1 hybrids

cd ~/sch_hae_scan

mkdir results/sim_f1_hybrid_reads

cd results/sim_f1_hybrid_reads

#get reference lists (minimal bovis or haem)
ln -s ../admixture/sh_ref.list
ln -s ../admixture/sb_ref.list


#get HAP1 and HAP2 fasta files
conda activate ~/sch_hae_scan/envs/bcftools
bgzip -c ../phasing/beagle.vcf >beagle.vcf.gz
tabix -p vcf beagle.vcf.gz

for INDIV in $(cat sh_ref.list sb_ref.list) ; do
  for HAP in 1 2; do

    echo $INDIV-h$HAP
    
    #generate haplotype fasta
    #conda activate ~/sch_hae_scan/envs/bcftools

    #bcftools consensus \
    #    -H $HAP \
    #    -M "?" \
    #    --sample $INDIV \
    #    -f ~/sch_hae_scan/data/schHae2_wmito.fa \
    #    beagle.vcf.gz \
    #    >$INDIV"_H"$HAP.fasta
    
    #conda deactivate


    # #simulate the reads
    # # conda activate ~/sch_hae_scan/envs/art

    # art_illumina \
    #     -sam \
    #     -i ~/sch_hae_scan/data/schHae2_wmito.fa \
    #     -p \
    #     -l 150 \
    #     -ss HS25 \
    #     -f 25 \
    #     -m 500 \
    #     -s 10 \
    #     -o $INDIV"_H"$HAP"_unsorted" &

    # # conda deactivate

    #samtools sort sam file
    #conda activate ~/sch_hae_scan/envs/samtools
    samtools sort -o $INDIV"_H"$HAP.bam --threads 16 $INDIV"_H"$HAP"_unsorted.sam"

    done
done 

#now randomly pair a sh with a sb
shuf sh_ref.list | head >sh_ref_shuf.list
shuf sb_ref.list | head >sb_ref_shuf.list
paste -d"," sh_ref_shuf.list sb_ref_shuf.list >paired.list


#and merge the same files for the hybrids (to be feed into gatk pipeline)

conda activate ~/sch_hae_scan/envs/samtools
I=0
for PAIR in $(cat paired.list); do
    SH=$(echo $PAIR | cut -f1 -d",")
    SB=$(echo $PAIR | cut -f2 -d",")

    OUT_BAM=in_silico_f1_n$I.bam

    samtools merge \
        --threads 16 \
        --output-fmt BAM \
        --reference ~/sch_hae_scan/data/schHae2_wmito.fa \
        $OUT_BAM
        $SH"_H"$(shuf -i 1-2 -n 1).bam
        $SB"_H"$(shuf -i 1-2 -n 1).bam

    I=$((I+1))

done
