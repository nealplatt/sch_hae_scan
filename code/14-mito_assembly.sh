#!/bin/bash

PROJ_DIR="/master/nplatt/sch_hae_scan"
RESULTS_DIR="$PROJ_DIR/results"
QSUB="qsub -V -cwd -S /bin/bash -q high_mem.q -j y -p -1023" 
CONDA="conda activate $PROJ_DIR/envs/mito_assembly"

conda activate $PROJ_DIR/envs/mito_assembly

mkdir $RESULTS_DIR/mito_assembly 

cp  $RESULTS_DIR/filter_genotypes/indvs_to_keep.list \
    $RESULTS_DIR/mito_assembly/samples.list
    
#just in case...add animal mt to database
get_organelle_config.py --add animal_mt

for SAMPLE in $(cat $RESULTS_DIR/mito_assembly/samples.list); do

    FQ1=$RESULTS_DIR/filtered_reads/"$SAMPLE"_filtered_R1.fq.gz 
    FQ2=$RESULTS_DIR/filtered_reads/"$SAMPLE"_filtered_R2.fq.gz
    FQX=$RESULTS_DIR/filtered_reads/"$SAMPLE"_filtered_RX.fq.gz

    CMD="get_organelle_from_reads.py \
        -1 $FQ1 \
        -2 $FQ2 \
        -u $FQX \
        -t 12 \
        -R 10 \
        -k 21,45,65,85,105 \
        -F animal_mt \
        -o $RESULTS_DIR/mito_assembly/$SAMPLE"
        
    echo "$CONDA; $CMD" | $QSUB -pe smp 12 -N mito_$SAMPLE -o results/logs/mito_$SAMPLE.log
done


