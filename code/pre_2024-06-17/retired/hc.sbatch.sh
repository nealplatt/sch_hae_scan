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


#sbatch --export=SAMPLE=$SAMPLE --job-name hc.$SAMPLE hc.sbatch.sh

# Launch serial code...

#get in proper dir and activate conda env
cd $SCRATCH/sch_hae_scan
# source $WORK/miniconda3/etc/profile.d/conda.sh
# conda init bash                      
# conda activate ascaris-01-filter_and_map

mkdir hc/$SAMPLE/

#run the large chromosomes in the background
for i in $(seq 1 8); do

    SCAFF=HiC_scaffold_"$i"

    ./gatk-4.2.0.0/gatk --java-options "-Xmx6g" HaplotypeCaller \
        --input mapped_reads/"$SAMPLE"_processed.bam \
        --output hc/$SAMPLE/"$SAMPLE"-"$SCAFF"_hc.vcf \
        --intervals $SCAFF \
        -reference data/ShV3_scaffolds.fa \
        --emit-ref-confidence GVCF \
        --native-pair-hmm-threads 4 \
        >logs/$SAMPLE"#"$SCAFF.log 2>&1 &
done

#run the smaller scaffolds in the foreground (serial)
for i in $(seq 9 563); do

    SCAFF=HiC_scaffold_"$i"

    ./gatk-4.2.0.0/gatk --java-options "-Xmx6g" HaplotypeCaller \
        --input mapped_reads/"$SAMPLE"_processed.bam \
        --output hc/$SAMPLE/"$SAMPLE"-"$SCAFF"_hc.vcf \
        --intervals $SCAFF \
        -reference data/ShV3_scaffolds.fa \
        --emit-ref-confidence GVCF \
        --native-pair-hmm-threads 4 \
        >logs/$SAMPLE"#"$SCAFF.log 2>&1

done

#do not exit until all major chroms are done
wait

#pause just a bit longer
sleep 1m