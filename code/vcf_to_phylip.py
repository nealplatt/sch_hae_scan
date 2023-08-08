import argparse
import vcf
import os
from collections import defaultdict
from multiprocessing import Pool, cpu_count
import time

#################################################
def alleles_to_iupac(h1, h2):
    """
    Returns the IUPAC code for two nucleotides.

    Parameters:
        h1 (str): The first nucleotide.
        h2 (str): The second nucleotide.

    Returns:
        The IUPAC code for the two nucleotides.

    """
    lookup = {"AG": "R", "GA": "R", "CT": "Y", "TC": "Y", "GC": "S", "CG": "S",
              "AT": "W", "TA": "W", "GT": "K", "TG": "K", "AC": "M", "CA": "M", 
              "AA": "A", "TT": "T", "CC": "C", "GG": "G", "..": "."}

    pair = h1 + h2
    return lookup.get(pair, "N")

#################################################
def vcf_to_phy(input_vcf, haploid, iupac, out_dir=None):

    # Open the VCF file
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))

    # Get the name of the PHYLIP file
    phylip_name = os.path.splitext(input_vcf)[0] + '.phy'

    # If the out_dir is provided, add it to the phylip_name and create the directory if it doesn't exist
    if out_dir:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        phylip_name = os.path.join(out_dir, os.path.basename(phylip_name))

    # Open the PHYLIP file
    haplotype_file = open(phylip_name, 'w')

    # Get the number of samples and the number of SNPs
    num_samples = len(vcf_reader.samples)
    num_snps = 0
    for record in vcf_reader:
        num_snps += 1

    # Get the reference and alternative alleles for each SNP
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))
    snps = []
    for record in vcf_reader:
        snps.append(record)

    # Get the sequence data for each sample into a single string (one per haplotype)
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))

    seq_data = defaultdict(lambda: defaultdict(str))

    for record in vcf_reader:
        for i in range(num_samples):
            sample = vcf_reader.samples[i]

            allele_1 = record.samples[i]['GT'][0]
            allele_2 = record.samples[i]['GT'][2]

            if allele_1 == str(0):
                allele_1 = str(record.REF)
            elif allele_1 == str(1):
                allele_1 = str(record.ALT[0])
            elif allele_1 == ".":
                allele_1 = "."

            if allele_2 == str(0):
                allele_2 = str(record.REF)
            elif allele_2 == str(1):
                allele_2 = str(record.ALT[0])
            elif allele_2 == ".":
                allele_2 = "."

            #If you are generating two haplotypes per sample
            if haploid:
                seq_data[sample]["1"] += allele_1
                seq_data[sample]["2"] += allele_2

            # If you are generating a single sequence per sample with heterozygous alleles as IUPAC
            elif iupac:
                allele = alleles_to_iupac(allele_1, allele_2)
                seq_data[sample]["IUPAC"] += allele


    if haploid:
        # Write the header for the PHYLIP file
        haplotype_file.write("{} {}\n".format(str(num_samples*2), num_snps))

        # Write the data the PHYLIP file
        for sample in vcf_reader.samples:
            haplotype_file.write("{}_1\t{}\n".format(sample, seq_data[sample]["1"]))
            haplotype_file.write("{}_2\t{}\n".format(sample, seq_data[sample]["2"]))    
    elif iupac:
        # Write the header for the PHYLIP file
        haplotype_file.write("{} {}\n".format(str(num_samples), num_snps))
        
        # Write the data the PHYLIP file
        for sample in vcf_reader.samples:
            haplotype_file.write("{}_1\t{}\n".format(sample, seq_data[sample]["IUPAC"]))

    # Close the PHYLIP file
    haplotype_file.close()


#################################################
# Create an argument parser
parser = argparse.ArgumentParser()

# Add the input VCF file argument
parser.add_argument('--input_vcf', help='Input a single VCF file. Ex. "locus_1.vcf"')

# Add the input VCF file list argument
parser.add_argument('--input_vcf_list', help='Input VCF file list. Ex. "loci_vcfs.list"')

# Add the output directory argument
parser.add_argument('--out_dir', help='Output directory for the PHYLIP file.  Defaults to the current directory.', default=os.getcwd(), )

# Add the number of threads
parser.add_argument('--threads', help='Number of threads. This only used when there are multiple VCF files provided with "--input_vcf_list" option', type=int, default=1)

# Require at least one of the new options to be used
group = parser.add_mutually_exclusive_group(required=True)
# Add the IUPAC option (default is False)

group.add_argument('--iupac', dest='iupac', action='store_true', help='Return the IUPAC code for each site. Ex. SampleA', default=False)

# Add the haploid option (default is False)
group.add_argument('--haploid', dest='haploid', action='store_true', help='Return a new sample for each haplotype. Ex. SampleA_1, SampleA_2', default=False)

# Parse the arguments
args = parser.parse_args()

#################################################
if args.input_vcf:
    # Call the vcf_to_phy function with the input VCF file and the output directory
    vcf_to_phy(args.input_vcf, args.haploid, args.iupac, args.out_dir)

elif args.input_vcf_list:
    # Open the list of VCF files
    with open(args.input_vcf_list, 'r') as vcf_list:
        # For each line in the file
        vcf_files = [line.strip() for line in vcf_list]

    # Determine the number of worker threads to use
    num_threads = min(args.threads, len(vcf_files))

    # Create a pool of worker processes and process the VCF files in parallel
    with Pool(num_threads) as p:
        # Use starmap_async to start the worker processes with the vcf_to_phy function and the list of VCF files as arguments
        results = p.starmap_async(vcf_to_phy, [(vcf_file, args.haploid, args.iupac, args.out_dir) for vcf_file in vcf_files], chunksize=1)

        # Check if the processing is complete, and print a message indicating the percentage completed and the number of VCF files remaining
        while not results.ready():
            # Calculate the total number of VCF files and the number of VCF files remaining to be processed
            num_files = len(vcf_files)
            remaining = results._number_left

            # Calculate the percentage of VCF files that have been processed
            percent_completed = float(((num_files-remaining)/num_files)*100)

            # Print a message indicating the percentage completed and the number of VCF files remaining
            print(f"{percent_completed:.2f}% completed. Processing {remaining} remaining VCF files...", end='\r')

            # Wait for half a second before checking again
            time.sleep(0.5)

        # Print a message indicating that the processing is complete
        time.sleep(2)
        print("Processing complete!")