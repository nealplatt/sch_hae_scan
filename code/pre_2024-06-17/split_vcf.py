import argparse
import vcf
import os
import tqdm


def split_vcf(input_file, num_sites, output_dir, output_csv):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create a list to store the data for the positions CSV file
    csvs = []

    # Create a defaultdict to store sites for each chromosome
    sites_per_chrom = {}

    # Read the VCF file
    vcf_reader = vcf.Reader(open(input_file, 'r'))

    # Process each record in the VCF file
    for record in tqdm.tqdm(vcf_reader):
        # Get the chromosome name
        chrom = record.CHROM

        # Add the record to the sites for the chromosome
        if chrom not in sites_per_chrom:
            sites_per_chrom[chrom] = []
        sites_per_chrom[chrom].append(record)

        # If the number of sites for the chromosome has reached the limit, write the output file
        if len(sites_per_chrom[chrom]) == num_sites:
            # Get the start and end positions of the split sites
            start_pos = sites_per_chrom[chrom][0].POS
            end_pos = sites_per_chrom[chrom][-1].POS

            # Create the output filename for the split VCF file
            out_filename = os.path.join(output_dir, "{}:{}-{}.vcf".format(chrom, start_pos, end_pos))

            # Write the split VCF file
            with open(out_filename, 'w') as outfile:
                vcf_writer = vcf.Writer(outfile, vcf_reader)
                for site in sites_per_chrom[chrom]:
                    vcf_writer.write_record(site)

            # Add the results to the list for the positions CSV file
            csvs.append((chrom, start_pos, end_pos, num_sites))

            # Clear the sites for the chromosome
            del sites_per_chrom[chrom]

    # Write the positions CSV file
    with open(output_csv, 'w') as outfile:
        outfile.write("#chrom,start,stop,num_sites\n")
        for result in csvs:
            outfile.write("{},{},{},{}\n".format(*result))

if __name__ == "__main__":
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Split a VCF file into smaller files')
    parser.add_argument('--input', '-i', type=str, required=True, help='Path to input VCF file')
    parser.add_argument('--num-sites', '-n', type=int, required=True, help='Number of sites per output VCF file')
    parser.add_argument('--output-dir', '-o', type=str, required=True, help='Output directory for output VCF file names')
    parser.add_argument('--output-csv', '-c', type=str, required=True, help='Name of the output CSV file')
    args = parser.parse_args()

    # Call the split_vcf function with the parsed arguments
    split_vcf(args.input, args.num_sites, args.output_dir, args.output_csv)

