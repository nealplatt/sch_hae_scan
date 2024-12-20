#!/usr/bin/bash

# Activate the conda environment for genome processing
conda activate scan-01-process_genome

# Set the project directory
PROJ_DIR=/master/nplatt/sch_hae_scan

# Create the project directory and subdirectories
mkdir -p ${PROJ_DIR}
cd ${PROJ_DIR}
mkdir results envs code bin data logs

# Prepare the genome
cd ${PROJ_DIR}/data

# Set the reference genome filename
REF="GCF_000699445.3_UoM_Shae.V3_genomic.fna"

# Download the reference genome dataset
datasets download genome accession GCF_000699445.3 --filename GCF_000699445.3.zip

# Unzip the downloaded genome dataset
unzip GCF_000699445.3.zip

# Create symbolic links to the reference genome and GFF files
ln -s ncbi_dataset/data/GCF_000699445.3/GCF_000699445.3_UoM_Shae.V3_genomic.fna .
ln -s ncbi_dataset/data/GCF_000699445.3/genomic.gff ./GCF_000699445.3_UoM_Shae.V3_genomic.gff

# Index the reference genome using BWA
bwa index $REF

# Activate the conda environment for filtering, mapping, and genotyping
conda activate scan-02-filter_map_genotype

# Create a sequence dictionary for the reference genome using GATK
$PROJ_DIR/bin/gatk-4.2.0.0/gatk CreateSequenceDictionary -R $REF

# Deactivate the conda environment
conda deactivate

# Index the reference genome using samtools
samtools faidx $REF > samtools.faidx.log

# Download SRA data
mkdir -p $PROJ_DIR/data/sra_data
cd $PROJ_DIR/data/sra_data

# Create an accession list of Schistosoma samples
cat <<EOL > SraAccList.txt
ERR103048 bovis_tanzania_ERR103048
ERR119612 guineensis_saotome_ERR119612
ERR119621 smargrebowiei_zambia_ERR119621
ERR119623 curassoni_senegal_ERR119623
ERR539851 bovis_unk_ERR539851
SRR433862 sh_egypt_SRR433862
ERR103051 matthei_zambia_ERR103051
ERR119613 intercalatum_drcongo_ERR119613
ERR119622 bovis_tanzania_ERR119622
ERR310940 margrebowiei_zambia_ERR310940
ERR539853 bovis_unk_ERR539853
SRR13579865 sh_egypt_SRR13579865
SRR13579866 sh_ivorycoast_SRR13579866
SRR13579867 sh_mali_SRR13579867
SRR13579868 sh_zanzibar_SRR13579868
SRR13579869 sh_scan_SRR13579869
SRR13579870 sh_guineabissau_SRR13579870
SRR13579871 sh_madagascar_SRR13579871
SRR13579872 sh_gambia_SRR13579872
SRR13579881 sh_senegal_SRR13579881
SRR13579882 sh_scan_SRR13579882
SRR13579883 sh_cameroon_SRR13579883
SRR8284792 sh_tzpem_SRR8284792
SRR8284793 sh_tzpem_SRR8284793
SRR8284794 sh_tzpem_SRR8284794
SRR8284791 sh_niger_SRR8284791
SRR8284786 sh_niger_SRR8284786
SRR8284797 sh_tzpem_SRR8284797
SRR8284796 sh_tzpem_SRR8284796
SRR8284795 sh_tzpem_SRR8284795
SRR8284788 sh_niger_SRR8284788
SRR8284787 sh_niger_SRR8284787
SRR8284790 sh_niger_SRR8284790
SRR8284789 sh_niger_SRR8284789
SRR13579873 bovis_spain_SRR13579873
SRR13579874 bovis_ethiopia_SRR13579874
SRR13579875 bovis_uganda_SRR13579875
SRR13579876 bovis_senegal_SRR13579876
SRR13579877 bovis_scan_SRR13579877
SRR13579878 bovis_keyna_SRR13579878
SRR7743803 f1sbsh_unk_SRR7743803
SRR7867226 bovis_tanzania_SRR7867226
SRR7867225 bovis_tanzania_SRR7867225
EOL

# Download and process each SRA sample
while read -r LINE; do
    # Extract SRA accession and sample name
    SRA=$(echo $LINE | cut -f1 -d" ")
    NAME=$(echo $LINE | cut -f2 -d" ")

    echo $NAME

    # Download the SRA file
    prefetch --max-size 400g --output-directory . $SRA

    # Convert SRA to FASTQ and compress
    fastq-dump --outdir ./$SRA --gzip --split-files ./$SRA/$SRA.sra &

    # Remove existing FASTQ files with the same name
    rm ~/sch_hae_scan/data/seq_data/"$NAME"_R*.fq.gz

    # Create symbolic links to the new FASTQ files
    ln -s ~/sch_hae_scan/data/sra_data/$SRA/"$SRA"_1.fastq.gz ~/sch_hae_scan/data/seq_data/"$NAME"_R1.fq.gz
    ln -s ~/sch_hae_scan/data/sra_data/$SRA/"$SRA"_2.fastq.gz ~/sch_hae_scan/data/seq_data/"$NAME"_R2.fq.gz

done < SraAccList.txt

##################################################################################
# Proceed to 02-filter_map_genotype.snake
##################################################################################
