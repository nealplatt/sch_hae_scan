conda activate scan_phylo

mkdir ~/sch_hae_scan/results/gene_trees
cd ~/sch_hae_scan/results/gene_trees


#remove outgroups and get major chrs (from maf05)
vcftools \
    --vcf ../phasing/beagle_maf05.vcf \
    --chr HiC_scaffold_1 \
    --chr HiC_scaffold_2 \
    --chr HiC_scaffold_3 \
    --chr HiC_scaffold_4 \
    --chr HiC_scaffold_5 \
    --chr HiC_scaffold_6 \
    --chr HiC_scaffold_7 \
    --chr HiC_scaffold_8 \
    --remove-indv  matthei_zambia_ERR103051 \
    --remove-indv  margrebowiei_zambia_ERR310940 \
    --remove-indv guineensis_saotome_ERR119612 \
    --remove-indv intercalatum_drcongo_ERR119613 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >beagle_maf05_no_outgroups_chrs.vcf

#get the outgroups and make sure the variant is fixed.
vcftools \
    --vcf ../filter_genotypes/maf05.vcf \
    --chr HiC_scaffold_1 \
    --chr HiC_scaffold_2 \
    --chr HiC_scaffold_3 \
    --chr HiC_scaffold_4 \
    --chr HiC_scaffold_5 \
    --chr HiC_scaffold_6 \
    --chr HiC_scaffold_7 \
    --chr HiC_scaffold_8 \
    --indv matthei_zambia_ERR103051 \
    --indv margrebowiei_zambia_ERR310940 \
    --max-maf 0 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >maf05_outgroups_fixed_chrs.vcf

#get only smarg/smatt gtd sites
vcftools \
    --vcf maf05_outgroups_fixed_chrs.vcf \
    --indv matthei_zambia_ERR103051 \
    --indv margrebowiei_zambia_ERR310940 \
    --missing-site \
    --stdout \
    >maf05_outgroups_fixed_chrs.missing-site

cat maf05_outgroups_fixed_chrs.missing-site  | awk '{if ($6<=0.5) print $1":"$2}' >sites_to_keep.list

#get same sites for each set of ingroup and outgroup samples
for VCF in beagle_maf05_no_outgroups_chrs.vcf maf05_outgroups_fixed_chrs.vcf; do
    vcftools \
        --vcf $VCF \
        --snps sites_to_keep.list \
        --recode \
        --recode-INFO-all \
        --stdout \
        >$(basename $VCF .vcf)"_outsites.vcf" &
done

#for whatever reason that doesn't work the best so I am doing another version to confirm its the same sites
# basically look at the same sites in each file via grep and counting occurences
grep -v "#" *_outsites.vcf | cut -f3 | sort | uniq -c | awk '{if ($1==2) print $2}' >confirm_sites_to_keep.list

for VCF in beagle_maf05_no_outgroups_chrs_outsites.vcf maf05_outgroups_fixed_chrs_outsites.vcf; do
    vcftools \
        --vcf $VCF \
        --snps confirm_sites_to_keep.list \
        --recode \
        --recode-INFO-all \
        --stdout \
        >$(basename $VCF .vcf)"_confirmed.vcf" &
done

#compress and index for merging
for VCF in $(ls *_confirmed.vcf); do
    echo $VCF
    #bgzip -c $VCF > $VCF.gz &
    tabix -p vcf $VCF.gz &
done

#now merge the vcf.gzs 
bcftools merge \
    -m snps \
    -o gene_tree_snps.vcf \
    -O v \
    beagle_maf05_no_outgroups_chrs_outsites_confirmed.vcf.gz \
    maf05_outgroups_fixed_chrs_outsites_confirmed.vcf.gz

bgzip -c gene_tree_snps.vcf >gene_tree_snps.vcf.gz
tabix -p vcf gene_tree_snps.vcf.gz


#split into equal windows of 500 snps.
mkdir snp_lists
for CHR in 1 2 3 4 5 6 7 8; do
     bcftools query -f'%CHROM\t%POS\n' --regions HiC_scaffold_${CHR} gene_tree_snps.vcf.gz >snp_lists/HiC_scaffold_${CHR}.list &
done

mkdir split_snp_lists
for CHR in 1 2 3 4 5 6 7 8; do
    echo $CHR
    split \
        --lines 500 \
        --suffix-length=6 \
        --additional-suffix=.list \
        --numeric-suffixes \
        snp_lists/HiC_scaffold_${CHR}.list \
        split_snp_lists/HiC_scaffold_${CHR}- &
done

#create vcf files
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -p -1023" 
mkdir logs
mkdir gene_vcfs

#split larger file into individual chroms to make it go much faster
for CHR in 1 2 3 4 5 6 7 8; do
    bcftools view gene_tree_snps.vcf.gz --regions HiC_scaffold_${CHR} -o HiC_scaffold_${CHR}.bcf -O b &
done

for CHR in 1 2 3 4 5 6 7 8; do
    SCAFF=HiC_scaffold_${CHR}
    for SNP_LIST in $(ls split_snp_lists/${SCAFF}-??????.list); do
        CMD="conda activate scan_phylo; bcftools view -T ${SNP_LIST} -O v ${SCAFF}.bcf >gene_vcfs/$(basename $SNP_LIST .list).vcf"
        
        echo $CMD  | $QSUB -pe smp 3 -N $(basename $SNP_LIST .list) -o logs/$(basename $SNP_LIST).log
    done
done

#haploidify each vcf file with custom script
mkdir gene_phys

for CHR in 1 2 3 4 5 6 7 8; do
    SCAFF=HiC_scaffold_${CHR}
    for VCF in $(ls gene_vcfs/${SCAFF}-??????.vcf); do
        #echo $VCF
        PHY=gene_phys/$(basename $VCF .vcf).phy
        python ./diploid_vcf_to_haploid_phy.py $VCF $PHY gene_tree_cords.csv
    done &
done

#create gene and bootstrap gene trees
mkdir gene_trees

for PHY in $(ls gene_phys/HiC_scaffold_?-??????.phy); do
    raxml-ng \
        --all \
        --msa $PHY \
        --msa-format PHYLIP \
        --model GTR+ASC_LEWIS \
        --tree pars{10} \
        --bs-trees 100 \
        --threads 10 \
        --prefix gene_trees/$(basename $PHY .phy) \
        --outgroup margrebowiei_zambia_ERR310940_h1,matthei_zambia_ERR103051_h1

done

#cat all gene trees

#remove low support branches
nw_ed  1KP-genetrees.tre 'i & b<=10' o > 1KP-genetrees-BS10.tre



###################
#run twist
###################
mkdir ~/sch_hae_scan/results/twisst
cd ~/sch_hae_scan/results/twisst

#run twist
cp ../gene_trees/genetrees.tre .
bcftools query -l ../gene_trees/gene_tree_snps.vcf >samples.list


 python ~/sch_hae_scan/bin/twisst/twisst.py
  -t genetrees.tre
  -outgroup 
  --groupsFile



###################
#svd quartets
###################

#convert vcf to phylip
cd ~/sch_hae_scan/bin
git clone https://github.com/edgardomortiz/vcf2phylip.git

cd ~/sch_hae_scan/results
mkdir svdq
cd svdq


#ld filtering
plink \
    --vcf ../gene_trees/gene_tree_snps.vcf \
    --allow-extra-chr \
    --double-id \
    --indep-pairwise 25 5 0.20 \
    --out ld_filtered

#Pruning complete.  6801863 of 7231736 variants removed.
vcftools \
    --vcf ../gene_trees/gene_tree_snps.vcf \
    --exclude ld_filtered.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >gene_tree_snps_filtered.vcf

python ../gene_trees/diploid_vcf_to_haploid_phy.py ../gene_trees/gene_tree_snps.vcf gene_tree_snps_hs.phy coords.txt

#convert phylip to nexus
run paup
toNEXUS fromfile=gene_tree_snps_filtered.phy format=RelPHYLIP tofile=gene_tree_snps_filtered.nex datatype=nucleotide interleaved=no replace=yes

exec gene_tree_snps_filtered.nex




begin sets;
    taxpartition groups =
        s_bovis_nigeria:             Sb_NG_au_1.2_h1 Sb_NG_au_1.2_h2 Sb_NG_au_2.10_h1 Sb_NG_au_2.10_h2 Sb_NG_au_2.13_h1 Sb_NG_au_2.13_h2 Sb_NG_au_2.5_h1 Sb_NG_au_2.5_h2 Sb_NG_au_2.6_h1 Sb_NG_au_2.6_h2 Sb_NG_be_1.10_h1 Sb_NG_be_1.10_h2 Sb_NG_be_1.3_h1 Sb_NG_be_1.3_h2 Sb_NG_be_1.5_h1 Sb_NG_be_1.5_h2 Sb_NG_en_1.1_h1 Sb_NG_en_1.1_h2,
        s_bovis_ethiopia:            bovis_ethiopia_SRR13579874_h1 bovis_ethiopia_SRR13579874_h2,
        s_bovis_kenya:               bovis_keyna_SRR13579878_h1 bovis_keyna_SRR13579878_h2,
        s_bovis_senegal:             bovis_senegal_SRR13579876_h1 bovis_senegal_SRR13579876_h2,
        s_bovis_tanzania:            bovis_tanzania_ERR103048_h1 bovis_tanzania_ERR103048_h2 bovis_tanzania_SRR7867225_h1 bovis_tanzania_SRR7867225_h2 bovis_tanzania_SRR7867226_h1 bovis_tanzania_SRR7867226_h2,
        s_bovis_cote_d_ivoire:       sbo_cdivoire_raviart_108_h1 sbo_cdivoire_raviart_108_h2 sbo_cdivoire_raviart_109_h1 sbo_cdivoire_raviart_109_h2 sbo_cdivoire_foro_106_h1 sbo_cdivoire_foro_106_h2 sbo_cdivoire_foro_116_h1 sbo_cdivoire_foro_116_h2 sbo_cdivoire_noumousso_107_h1 sbo_cdivoire_noumousso_107_h2 sbo_cdivoire_noumousso_117_h1 sbo_cdivoire_noumousso_117_h2 ssp_cdivoire_foro_101_h1 ssp_cdivoire_foro_101_h2 ssp_cdivoire_foro_103_h1 ssp_cdivoire_foro_103_h2 ssp_cdivoire_noumousso_88_h1 ssp_cdivoire_noumousso_88_h2,
        s_bovis_uganda:              sbo_uganda_apac13_35_h1 sbo_uganda_apac13_35_h2,
        s_bovis_niger:               sbo_niger_libore_154_h1 sbo_niger_libore_154_h2 sbo_niger_libore_155_h1 sbo_niger_libore_155_h2 sbo_niger_libore_160_h1 sbo_niger_libore_160_h2 sbo_niger_libore_162_h1 sbo_niger_libore_162_h2 sbo_niger_libore_163_h1 sbo_niger_libore_163_h2 sbo_niger_libore_164_h1 sbo_niger_libore_164_h2 sbo_niger_libore_170_h1 sbo_niger_libore_170_h2 sbo_niger_libore_171_h1 sbo_niger_libore_171_h2 sbo_niger_libore_172_h1 sbo_niger_libore_172_h2 sbo_niger_libore_175_h1 sbo_niger_libore_175_h2 sbo_niger_libore_176_h1 sbo_niger_libore_176_h2,
        s_haematobium_nigeria:       Sh_NG_eb_6_2_h1 Sh_NG_eb_6_2_h2 Sh_NG_ed_1_3_h1 Sh_NG_ed_1_3_h2 Sh_NG_ed_3_2_h1 Sh_NG_ed_3_2_h2 Sh_NG_kb_2_1_h1 Sh_NG_kb_2_1_h2 Sh_NG_kw_1_10_h1 Sh_NG_kw_1_10_h2 Sh_NG_os_1_4_h1 Sh_NG_os_1_4_h2 sh_nigeria_SRR13579869_h1 sh_nigeria_SRR13579869_h2,
        s_haematobium_cote_d_ivoire: sh_ivorycoast_SRR13579866_h1 sh_ivorycoast_SRR13579866_h2 sha_cdivoire_kongobo_115_h1 sha_cdivoire_kongobo_115_h2 ssp_cdivoire_allokokro_102_h1 ssp_cdivoire_allokokro_102_h2 ssp_cdivoire_allokokro_104_h1 ssp_cdivoire_allokokro_104_h2 ssp_cdivoire_kongobo_95_h1 ssp_cdivoire_kongobo_95_h2 ssp_cdivoire_kongobo_96_h1 ssp_cdivoire_kongobo_96_h2 ssp_cdivoire_kongobo_97_h1 ssp_cdivoire_kongobo_97_h2 ssp_cdivoire_kongobo_98_h1 ssp_cdivoire_kongobo_98_h2 ssp_cdivoire_linguebo_100_h1 ssp_cdivoire_linguebo_100_h2 ssp_cdivoire_linguebo_105_h1 ssp_cdivoire_linguebo_105_h2 ssp_cdivoire_noumousso_94_h1 ssp_cdivoire_noumousso_94_h2 ssp_cdivoire_raviart_110_h1 ssp_cdivoire_raviart_110_h2 ssp_cdivoire_raviart_111_h1 ssp_cdivoire_raviart_111_h2 ssp_cdivoire_raviart_113_h1 ssp_cdivoire_raviart_113_h2 ssp_cdivoire_raviart_114_h1 ssp_cdivoire_raviart_114_h2 ssp_cdivoire_raviart_89_h1 ssp_cdivoire_raviart_89_h2 ssp_cdivoire_raviart_90_h1 ssp_cdivoire_raviart_90_h2 ssp_cdivoire_raviart_91_h1 ssp_cdivoire_raviart_91_h2 ssp_cdivoire_raviart_92_h1 ssp_cdivoire_raviart_92_h2 ssp_cdivoire_raviart_99_h1 ssp_cdivoire_raviart_99_h2,
        s_haematobium_niger:         sh_niger_SRR8284786_h1 sh_niger_SRR8284786_h2 sh_niger_SRR8284787_h1 sh_niger_SRR8284787_h2 sh_niger_SRR8284788_h1 sh_niger_SRR8284788_h2 sh_niger_SRR8284789_h1 sh_niger_SRR8284789_h2 sh_niger_SRR8284790_h1 sh_niger_SRR8284790_h2 sh_niger_SRR8284791_h1 sh_niger_SRR8284791_h2 sha_niger_libore_138_h1 sha_niger_libore_138_h2 sha_niger_libore_139_h1 sha_niger_libore_139_h2 sha_niger_libore_140_h1 sha_niger_libore_140_h2 sha_niger_libore_141_h1 sha_niger_libore_141_h2 sha_niger_libore_142_h1 sha_niger_libore_142_h2 sha_niger_libore_143_h1 sha_niger_libore_143_h2 sha_niger_libore_146_h1 sha_niger_libore_146_h2 sha_niger_libore_147_h1 sha_niger_libore_147_h2 sha_niger_libore_148_h1 sha_niger_libore_148_h2 sha_niger_libore_149_h1 sha_niger_libore_149_h2 sha_niger_libore_150_h1 sha_niger_libore_150_h2 sha_niger_libore_151_h1 sha_niger_libore_151_h2 ssp_niger_libore_156_h1 ssp_niger_libore_156_h2 ssp_niger_libore_157_h1 ssp_niger_libore_157_h2 ssp_niger_libore_159_h1 ssp_niger_libore_159_h2 ssp_niger_libore_165_h1 ssp_niger_libore_165_h2 ssp_niger_libore_166_h1 ssp_niger_libore_166_h2 ssp_niger_libore_167_h1 ssp_niger_libore_167_h2 ssp_niger_libore_168_h1 ssp_niger_libore_168_h2 ssp_niger_libore_169_h1 ssp_niger_libore_169_h2 ssp_niger_libore_173_h1 ssp_niger_libore_173_h2,
        s_haematobium_uganda:        sbo_uganda_runga_44_h1 sbo_uganda_runga_44_h2,
        s_haematobium_cameroon:      sh_cameroon_SRR13579883_h1 sh_cameroon_SRR13579883_h2,
        s_haematobium_egypt:         sh_egypt_SRR13579865_h1 sh_egypt_SRR13579865_h2,
        s_haematobium_gambia:        sh_gambia_SRR13579872_h1 sh_gambia_SRR13579872_h2,
        s_haematobium_guineabissau:  sh_guineabissau_SRR13579870_h1 sh_guineabissau_SRR13579870_h2,
        s_haematobium_mali:          sh_mali_SRR13579867_h1 sh_mali_SRR13579867_h2,
        s_haematobium_senegal:       sh_senegal_SRR13579881_h1 sh_senegal_SRR13579881_h2,
        s_haematobium_zanzibar:      sh_tzpem_SRR8284792_h1 sh_tzpem_SRR8284792_h2 sh_tzpem_SRR8284794_h1 sh_tzpem_SRR8284794_h2 sh_tzpem_SRR8284796_h1 sh_tzpem_SRR8284796_h2 sh_tzpem_SRR8284797_h1 sh_tzpem_SRR8284797_h2 sh_zanzibar_SRR13579868_h1 sh_zanzibar_SRR13579868_h2 sha_pemba_uwandani_1_h1 sha_pemba_uwandani_1_h2 sha_pemba_uwandani_10_h1 sha_pemba_uwandani_10_h2 sha_pemba_uwandani_11_h1 sha_pemba_uwandani_11_h2 sha_pemba_uwandani_12_h1 sha_pemba_uwandani_12_h2 sha_pemba_uwandani_13_h1 sha_pemba_uwandani_13_h2 sha_pemba_uwandani_15_h1 sha_pemba_uwandani_15_h2 sha_pemba_uwandani_2_h1 sha_pemba_uwandani_2_h2 sha_pemba_uwandani_5_h1 sha_pemba_uwandani_5_h2 sha_pemba_uwandani_7_h1 sha_pemba_uwandani_7_h2 sha_unguja_kinyasini_16_h1 sha_unguja_kinyasini_16_h2 sha_unguja_kinyasini_19_h1 sha_unguja_kinyasini_19_h2 sha_unguja_kinyasini_23_h1 sha_unguja_kinyasini_23_h2 sha_unguja_kinyasini_25_h1 sha_unguja_kinyasini_25_h2 sha_unguja_kinyasini_26_h1 sha_unguja_kinyasini_26_h2 sha_unguja_kinyasini_27_h1 sha_unguja_kinyasini_27_h2 sha_unguja_kinyasini_28_h1 sha_unguja_kinyasini_28_h2 sha_unguja_kinyasini_29_h1 sha_unguja_kinyasini_29_h2 sha_unguja_kinyasini_30_h1 sha_unguja_kinyasini_30_h2,
        s_haematobium_angola:        sha_angola_cota_57_h1 sha_angola_cota_57_h2 sha_angola_cota_58_h1 sha_angola_cota_58_h2 sha_angola_cota_59_h1 sha_angola_cota_59_h2 sha_angola_icau_60_h1 sha_angola_icau_60_h2 sha_angola_icau_61_h1 sha_angola_icau_61_h2 sha_angola_icau_62_h1 sha_angola_icau_62_h2 sha_angola_icau_63_h1 sha_angola_icau_63_h2 sha_angola_icau_64_h1 sha_angola_icau_64_h2,
        s_haematobium_kenya:         sha_kenya_unk_49_h1 sha_kenya_unk_49_h2,
        s_haematobium_liberia:       sha_liberia_dormeyan_54_h1 sha_liberia_dormeyan_54_h2 sha_liberia_joekpenmue_55_h1 sha_liberia_joekpenmue_55_h2,
        s_haematobium_madagascar:    sha_madag_belesalampy_73_h1 sha_madag_belesalampy_73_h2 sha_madag_belesalampy_74_h1 sha_madag_belesalampy_74_h2 sha_madag_belesalampy_75_h1 sha_madag_belesalampy_75_h2 sha_madag_belesalampy_76_h1 sha_madag_belesalampy_76_h2 sha_madag_belesalampy_77_h1 sha_madag_belesalampy_77_h2 sha_madag_belesalampy_78_h1 sha_madag_belesalampy_78_h2 sha_madag_belesalampy_80_h1 sha_madag_belesalampy_80_h2 sha_madag_belesalampy_81_h1 sha_madag_belesalampy_81_h2 sha_madag_belesalampy_82_h1 sha_madag_belesalampy_82_h2 sha_madag_belesalampy_83_h1 sha_madag_belesalampy_83_h2 sha_madag_belesalampy_84_h1 sha_madag_belesalampy_84_h2 sha_madag_belesalampy_87_h1 sha_madag_belesalampy_87_h2,
        s_haematobium_namibia:       sha_namib_mayenzere_134_h1 sha_namib_mayenzere_134_h2,
        s_haematobium_sudan:         sha_sudan_canal4_125_h1 sha_sudan_canal4_125_h2 sha_sudan_canal4_126_h1 sha_sudan_canal4_126_h2 sha_sudan_canal4_127_h1 sha_sudan_canal4_127_h2 sha_sudan_canal4_132_h1 sha_sudan_canal4_132_h2 sha_sudan_schoole_120_h1 sha_sudan_schoole_120_h2 sha_sudan_schoolh_123_h1 sha_sudan_schoolh_123_h2 sha_sudan_schooli_129_h1 sha_sudan_schooli_129_h2 sha_sudan_schooli_130_h1 sha_sudan_schooli_130_h2,
        s_haematobium_swaziland:     sha_swaz_mkhuzweni_177_h1 sha_swaz_mkhuzweni_177_h2 sha_swaz_mkhuzweni_178_h1 sha_swaz_mkhuzweni_178_h2 sha_swaz_mkhuzweni_179_h1 sha_swaz_mkhuzweni_179_h2 sha_swaz_mkhuzweni_180_h1 sha_swaz_mkhuzweni_180_h2 sha_swaz_mkhuzweni_182_h1 sha_swaz_mkhuzweni_182_h2 sha_swaz_mkhuzweni_184_h1 sha_swaz_mkhuzweni_184_h2 sha_swaz_njojane_188_h1 sha_swaz_njojane_188_h2 sha_swaz_njojane_190_h1 sha_swaz_njojane_190_h2 sha_swaz_njojane_191_h1 sha_swaz_njojane_191_h2 sha_swaz_qomintaba_185_h1 sha_swaz_qomintaba_185_h2 sha_swaz_unk_186_h1 sha_swaz_unk_186_h2 sha_swaz_unk_187_h1 sha_swaz_unk_187_h2,
        s_haematobium_zambia:        sha_zambia_lishiko_66_h1 sha_zambia_lishiko_66_h2 ssp_zambia_kafue_71_h1 ssp_zambia_kafue_71_h2,
        s_margrebowieie_na:          margrebowiei_zambia_ERR310940_h1,
        s_matthei_zambia:            matthei_zambia_ERR103051_h1;
end; 

exec taxpartition.nex


svdq evalq=random nquartets=1 000 000 taxpartition=groups bootstrap=standard nthreads=64;
svdq evalq=random nquartets=1000000 taxpartition=groups bootstrap=standard nreps=100 nthreads=192 showScores=yes;

svdq evalq=random nquartets=1000000 taxpartition=groups nthreads=192 showScores=yes bootstrap=standard nreps=100;


begin paup;
    outgroup snakespecies.Agkistrodon;
    set outroot=monophyl;
    taxcolor green      snakespecies.S.c._catenatus;
    [taxcolor (1 .7 1)      snakespecies.S.c._catenatus;]
    taxcolor dkgreen    snakespecies.S.m._miliarius;
    taxcolor turquoise  snakespecies.S.m._barbouri;
    taxcolor orange     snakespecies.S.m._streckeri;
    taxcolor red        snakespecies.S.c._edwardsii;
    taxcolor blue       snakespecies.S.c._tergeminus;
    taxcolor black      snakespecies.Agkistrodon;
end;

svdq evalq=random nquartets=1000 taxpartition=groups nthreads=64 treefile=svdq.trees;
savetrees file=test.tre

svdq evalq=random nquartets=1000 taxpartition=groups bootstrap=standard nreps=10 nthreads=64 treefile=svdq.trees;
savetrees file=boot.tre





svdq evalq=random nquartets=1000000 taxpartition=groups bootstrap=standard nreps=100 nthreads=64 treefile=svdq.trees;
savetrees file=svdq_boot.trees


nw_ed  HiC_scaffold_1-000008.raxml.support 'i & b<=50' o > test.tre


cat *.support >support.trees
nw_ed  support.trees 'i & b<=10' o >10.trees; nw_ed  support.trees 'i & b<=25' o >25.trees;



mkdir beast
cd beast

#randomly sample genes
ls ../gene_trees/split_snp_lists/*.list | shuf | head -n 500 >random_loci.list

mkdir random_loci
for LOCUS in $(cat random_loci.list); do
    cp ${LOCUS} random_loci/
done

mkdir svdq_validation
cat random_loci/*.list >svdq_validation/snps.list
#create vcf of only these snps for SVDQ validation

cd svdq_validation
cp ../../svdq/taxpartition.nex .

#create phylip file from vcf containing only target sites
vcftools \
    --vcf ../../gene_trees/gene_tree_snps.vcf \
    --positions snps.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >reduced_snps.vcf

conda activate diploid_vcf_to_haploid_vcf
python ../../gene_trees/diploid_vcf_to_haploid_phy.py \
    reduced_snps.vcf \
    reduced_snps_hs.phy \
    reduced_snps_coords.txt
conda deactivate

~/sch_hae_scan/bin/paup4a168_centos64
toNEXUS fromfile=reduced_snps_hs.phy format=RelPHYLIP tofile=reduced_snps_hs.phy.nex datatype=nucleotide interleaved=no replace=yes

exec reduced_snps_hs.phy.nex
exec taxpartition.nex

svdq evalq=random nquartets=1000000 taxpartition=groups bootstrap=standard nreps=100 nthreads=64 treefile=reduced_svdq.trees;
savetrees file=svdq_boot.trees





#cp snp file
cp ../gene_trees/gene_tree_snps.vcf .

#convert to phylip
~/sch_hae_scan/bin/vcf2phylip/vcf2phylip.py \
    -i gene_tree_snps_f

#run paup
~/sch_hae_scan/bin/paup4a168_centos64

exec gene_tree_snps.min4.nexus
svdq evalq=random nquartets=1000000 taxpartition=none bootstrap=standard nreps=100 nthreads=192 showScores=no speciesTree=no
savetrees file=svdq_boot.trees


nw_ed  support.trees 'i & b<=50' o >50.trees; nw_ed  support.trees 'i & b<=25' o >25.trees;

cat pca_df.csv | cut -f1,27,37  -d"," >../svdq/sample_info.csv


#################################

conda activate scan_phylo

mkdir ~/sch_hae_scan/results/svdq
cd ~/sch_hae_scan/results/svdq


vcftools \
    --vcf ../phasing/beagle_maf05.vcf \
    --chr HiC_scaffold_1 \
    --chr HiC_scaffold_2 \
    --chr HiC_scaffold_3 \
    --chr HiC_scaffold_4 \
    --chr HiC_scaffold_5 \
    --chr HiC_scaffold_6 \
    --chr HiC_scaffold_7 \
    --chr HiC_scaffold_8 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >beagle_maf05_chrs.vcf


#ld filtering
plink \
    --vcf ../phasing/beagle_maf05.vcf \
    --allow-extra-chr \
    --double-id \
    --indep-pairwise 25 5 0.20 \
    --out ld_filtered
#Pruning complete.  7230165 of 7718124 variants removed.

vcftools \
    --vcf ../phasing/beagle_maf05.vcf \
    --exclude ld_filtered.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >beagle_maf05_ld_filtered.vcf

vcftools \
    --vcf beagle_maf05_ld_filtered.vcf \
    --chr HiC_scaffold_1 \
    --chr HiC_scaffold_2 \
    --chr HiC_scaffold_3 \
    --chr HiC_scaffold_4 \
    --chr HiC_scaffold_5 \
    --chr HiC_scaffold_6 \
    --chr HiC_scaffold_7 \
    --chr HiC_scaffold_8 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >beagle_maf05_ld_filtered_chrs.vcf

#convert vcf to phylip
../../bin/vcf2phylip/vcf2phylip.py --input beagle_maf05_ld_filtered_chrs.vcf

#convert phylip to nexus
#run paup
~/sch_hae_scan/bin/paup4a168_centos64
toNEXUS fromfile=beagle_maf05_ld_filtered_chrs.min4.phy format=RelPHYLIP tofile=svdq.nex datatype=nucleotide interleaved=no replace=yes

exec svdq.nex


exec svdq.nex
svdq evalq=random nquartets=1000000 taxpartition=none bootstrap=standard nreps=100 nthreads=192 showScores=no speciesTree=no
savetrees file=svdq_boot.trees


