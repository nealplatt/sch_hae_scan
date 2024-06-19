conda activate mappable_reads

cd /master/nplatt/sch_hae_scan/results/mappable_regions

ln -s ../../data/GCF_000699445.3_UoM_Shae.V3_genomic.fna sh_v3.fa
ln -s ../../data/GCA_944470425.1_tdSchBovi1.1_genomic.fna sb.fa

art_illumina \
    --errfree \
    --fcov 100 \
    --out sh_100x \
    --paired \
    --in sh_v3.fa \
    --len 150 \
    --mflen 500 \
    --sdev 50 \
    >sh_art.log 2>&1 &

art_illumina \
    --errfree \
    --fcov 100 \
    --out sb_100x \
    --paired \
    --in sb.fa \
    --len 150 \
    --mflen 500 \
    --sdev 50 \
    >sb_art.log 2>&1 &

#map sh and sb reads back to sh_v3 genome (parameters?)
bbmap.sh \
    -ref=sh_v3.fa \
    -in=sb_100x1.fq \
    -in2=sb_100x2.fq \
    vslow \
    -threads=192 \
    ambig=toss \
    interleaved=false \
    -Xmx24g \
    -eoom \
    -out=sb_vs_shv3.bam \
    minid=0.8

bbmap.sh \
    -ref=sh_v3.fa \
    -in=sh_100x1.fq \
    -in2=sh_100x2.fq \
    vslow \
    -threads=192 \
    ambig=toss \
    interleaved=false \
    -Xmx24g \
    -eoom \
    -out=sh_vs_shv3.bam \
    minid=0.8

bbmap.sh \
    -ref=sb.fa \
    -in=sb_100x1.fq \
    -in2=sb_100x2.fq \
    vslow \
    -threads=192 \
    ambig=toss \
    interleaved=false \
    -Xmx24g \
    -eoom \
    -out=sb_vs_sb.bam \
    minid=0.8

bbmap.sh \
    -ref=sb.fa \
    -in=sh_100x1.fq \
    -in2=sh_100x2.fq \
    vslow \
    -threads=192 \
    ambig=toss \
    interleaved=false \
    -Xmx24g \
    -eoom \
    -out=sh_vs_sb.bam \
    minid=0.8

#run mosdepth to get areas with greater than 50x coverage
for COMBO in sb_vs_sb sh_vs_sb; do
    echo ${COMBO}
    samtools sort --threads 12 -o ${COMBO}_sorted.bam ${COMBO}.bam
    samtools index ${COMBO}_sorted.bam

    mosdepth -t 4 ${COMBO} ${COMBO}_sorted.bam
    zcat ${COMBO}.per-base.bed.gz } | awk 'BEGIN {FS="\t"}; {if ($4>50) print $0}' | awk 'BEGIN {FS="\t"}; {if ($4<150) print $0}' >${COMBO}_gt50x_lt150x.bed
    bedtools merge -i ${COMBO}_gt50x_lt150x.bed >${COMBO}_gt50x_lt150x_merged.bed
done


#find overlab between sh and sb intervals
#these are the white-list regions
bedtools intersect -a sh_vs_shv3_gt50x_lt150x_merged.bed -b sb_vs_sh_gt50x_lt150x_merged.bed >shv3_mappable_regions.bed
bedtools intersect -a sh_vs_sb_gt50x_lt150x_merged.bed -b sb_vs_sb_gt50x_lt150x_merged.bed >sb_mappable_regions.bed


# cat sh_vs_shv3_gt50x_lt150x_merged.bed | awk 'BEGIN {FS="\t"}; {SUM += $3-$2} END {print SUM}'
# 363840545 (89.8%)
#  cat sb_vs_sh_gt50x_lt150x_merged.bed | awk 'BEGIN {FS="\t"}; {SUM += $3-$2} END {print SUM}'
# 297534874 (73.4%)
# cat mappable_regions.bed | awk 'BEGIN {FS="\t"}; {SUM += $3-$2} END {print SUM}'
# 296509226 (73.2%)