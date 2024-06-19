import vcf
from tqdm import tqdm
import sys
import getopt


def count_vcf_records(vcf_file):
    count = 0
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                count += 1
    return count


def read_populations(population_file):
    pop_map = {}
    with open(population_file, 'r') as f:
        for line in f:
            sample, population = line.strip().split('\t')
            if population not in pop_map:
                pop_map[population] = []
            pop_map[population].append(sample)
    return pop_map

def process_vcf(vcf_file, pop_map, output_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    
    # Create a dictionary mapping samples to populations
    sample_to_pop = {ind: pop for pop, inds in pop_map.items() for ind in inds}
    
    with open(output_file, 'w') as out:
        # Write population names as header
        out.write(' '.join(pop_map.keys()) + '\n')
        
        total_iterations = count_vcf_records(vcf_file)

        # Process each record in the VCF file
        for record in tqdm(vcf_reader, total=total_iterations):
            allele_counts = {pop: [0, 0] for pop in pop_map}
            
            for sample in record.samples:
                if sample.called:
                    pop = sample_to_pop[sample.sample]
                    ref_count = sample['GT'].count('0')
                    alt_count = sample['GT'].count('1')
                    allele_counts[pop][0] += ref_count
                    allele_counts[pop][1] += alt_count
            
            # Write SNV information to the output file
            line_counts = [f"{counts[0]},{counts[1]}" for counts in allele_counts.values()]
            out.write(' '.join(line_counts) + '\n')

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "v:p:o:", ["vcf_file=", "pop_map=", "output_file="])
    except getopt.GetoptError:
        print("Usage: vcf2treemix.py -v <vcf_file> -p <pop_map> -o <output_file>")
        sys.exit(2)

    vcf_file = None
    pop_map = None
    output_file = None

    for opt, arg in opts:
        if opt in ("-v", "--vcf_file"):
            vcf_file = arg
        elif opt in ("-p", "--pop_map"):
            pop_map = arg
        elif opt in ("-o", "--output_file"):
            output_file = arg

    if vcf_file is None or pop_map is None or output_file is None:
        print("Usage: vcf2treemix.py -v <vcf_file> -p <pop_map> -o <output_file>")
        sys.exit(2)

    pop_map_data = read_populations(pop_map)
    process_vcf(vcf_file, pop_map_data, output_file)

if __name__ == "__main__":
    main(sys.argv[1:])
