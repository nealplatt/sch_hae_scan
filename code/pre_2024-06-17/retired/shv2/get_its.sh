cd /master/nplatt/sch_hae_scan

#mkdir results/its


REF_FAS="/master/nplatt/sch_hae_scan/data/schHae2_wmito.fa"
RESULTS_DIR="/master/nplatt/sch_hae_scan/results"

#Cattle as natural host for Schistosoma haematobium (Bilharz, 1852) Weinland, \
#   1858 x Schistosoma bovis Sonsino, 1876 interactions, with new cercarial \
#   emergence and genetic patterns

##MT158872.1 vs shv2

#cattle profile MT158872

#NW_023366102.1:1508664-1509323

#get its region
conda activate envs/vcftools

vcftools \
	--vcf $RESULTS_DIR/filter_genotypes/annotated_snps.vcf \
	--chr NW_023366102.1 \
	--from-bp 1508664 \
	--to-bp 1509323 \
	--recode \
	--recode-INFO-all \
	--stdout \
	>$RESULTS_DIR/its/its.vcf

#get the vcf for this contig only
#get the fasta file for this contig only
samtools faidx

cat ref.fa | vcf-consensus your_file.vcf.gz > out.fa

#physically phase?

#extract fasta


#phylogeny


#categorize as hh/hb/bb
