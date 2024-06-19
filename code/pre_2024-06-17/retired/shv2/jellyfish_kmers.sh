#kmer species id

cd ~/sch_hae_scan

mkdir results/kmers

cd results/kmers


#read in filtered R1s for each sample
conda activate ~/sch_hae_scan/envs/jellyfish

for R1 in $(ls ../filtered_reads/*filtered_R1.fq.gz); do

    SAMPLE=$(basename $R1 _filtered_R1.fq.gz)
    echo $SAMPLE

    jellyfish count \
        --mer-len 21 \
        --size 100M \
        --threads 16 \
        --canonical \
        --output $SAMPLE"_21kmer_counts.jf" \
        <(zcat $R1)

done

# #(add in silicos as well)
# for R1_H1 in $(ls ../sim_f1_hybrid_reads/*_H1_unsorted1.fq); do

#     SAMPLE=$(basename R1 _H1_unsorted1.fq)
#     R1_H2=$(echo $R1_H1 | sed 's/H1/H2/g')

#     jellyfish count \
#         --mer-len 21 \
#         --size 100M \
#         --threads 16 \
#         --canonical \
#         --output $SAMPLE"_21kmer_counts.jf \
#         <(zcat $R1_H1) <(zcat $R1_H2)

# done


#
jellyfish dump mer_counts.jf > mer_counts_dumps.fa
