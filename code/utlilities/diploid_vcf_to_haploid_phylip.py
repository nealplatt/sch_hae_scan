import os
import vcf
import re
from collections import defaultdict
import sys
from tqdm import tqdm
from pathlib import Path
import glob

vcf_file_list=sys.argv[1]
o_cord_file=sys.argv[2]

vcfs = glob.glob(i_vcf_regex)

for i_vcf_file in tqdm(vcfs):
    o_phy_file = "{}.phylip".format(i_vcf_file)

    vcf_reader = vcf.Reader(filename=i_vcf_file)
    samples = vcf_reader.samples

    # Defining a dict
    seqs=defaultdict(str)

    chrom="na"
    start=1e12
    stop=-1

    for record in vcf_reader:
        pos=record.POS
        if record.POS < start:
            start = record.POS

        if record.POS > stop:
            stop = record.POS

        chrom = record.CHROM

        for sample in samples:
            ref_allele=str(record.REF)
            alt_allele=str(record.ALT[0])

            h1 = re.split("\\|", record.genotype(sample)['GT'].replace("/", "|"))[0]
            h2 = re.split("\\|", record.genotype(sample)['GT'].replace("/", "|"))[1]

            if h1 == "0":
                h1=ref_allele
            if h1 == "1":
                h1=alt_allele

            if h2 == "0":
                h2=ref_allele
            if h2 == "1":
                h2=alt_allele

            if h1 == ".":
                h1="?"
            if h2 == ".":
                h2 = "?"

            seqs["{}_h1".format(sample)]+=h1
            seqs["{}_h2".format(sample)]+=h2
            # seqs["{}_h1".format(sample)].append(h1)
            # seqs["{}_h2".format(sample)].append(h2)

    num_samples = len(seqs.keys())
    seq_len = len(seqs["{}_h1".format(samples[0])])


    num_samples = len(seqs.keys())
    seq_len = len(seqs["{}_h1".format(samples[0])])


    with open(o_phy_file, 'w') as outf:
        outf.write("{} {}\n".format(num_samples, seq_len))
        for h_id in seqs:
            seq=seqs[h_id]
            outf.write("{}\t{}\n".format(h_id, seq))

    with open(o_cord_file, 'a') as outf:
        prefix = os.path.basename(i_vcf_file).replace(".vcf", "")
        outf.write("{},{},{},{}\n".format(prefix, chrom, start, stop))