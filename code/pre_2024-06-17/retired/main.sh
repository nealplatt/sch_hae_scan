#!/usr/bin/bash
conda activate scan-00-main

cd sch_hae_scan

mkdir -p data/seq_data
mkdir results code sandbox envs

#set up 
# - install gatk
# - git clone
# - make config file
#

#get sequence data from admera into data/seq_data

#get sra data
code/get_sra_data.sh

#get haematobium assembly
code/01-get_shv3_genome_and_sra_data.sh

#filter, map, and genotype with snakemake
snakemake \
  --printshellcmds \
  --use-conda \
  --cluster 'qsub -V -cwd -S /bin/bash -pe smp {threads} -o {log}.log -j y -q all.q' \
  --jobs 100 \
  --latency-wait 200 \
  --keep-going \
  --snakefile code/02-filter_map_genotype.snake \
  --rerun-incomplete 

#filter genotypes with snakemake
code/filter_gvcf.sh
    
#phase the variants
code/phasing.sh

################################################################################
#start running in jupyter notebooks
################################################################################
#start jupyter lab session
conda activate envs/jupyter_postproc

#this notebook assembles the mitochondrial genomes for all samples in the
#  filtered dataset, it also genotypes the ITS fragments often used in 
#  single gene PCR identification of hybrids.  The goals is to see how many and
#  which of our samples would be considered "hybrids" using traditional methods
code/notebooks/mito_and_its.ipynb

#this notebook looks at ancestry across the genome using multiple methods
#  including loter (ancestry assignment) and stats like D, F3, etc.  Some of
#  these stats are cacluated from the entire genome, and others are examined
#  in a slinding window *across* the genome
code/notebooks/ancestry.ipynb
