import os
import glob
from pathlib import Path

proj_dir="/master/nplatt/sch_hae_scan"
twisst_dir = f"{proj_dir}/results/twisst"
os.chdir(twisst_dir)

# Get a list of input files
input_files = glob.glob("gene_phylips/NC_067*.phy")
samples = [os.path.basename(f).replace('.phy', '') for f in input_files]

localrules: 
    all, 
    
# Run the workflow on all input files
rule all:
    input:
        expand(twisst_dir + "/gene_trees/{sample}.raxml.{ext}", sample=samples, ext=["bestTree"]),
        #expand(twisst_dir + "/gene_trees/{sample}.raxml.{ext}", sample=samples, ext=["bestTree", "support"]),
        f"{twisst_dir}/gene_trees.nwk",
        #twisst_dir + "/gene_trees_BS10.nwk"

# # Rule to run RAxML-NG on input files
# rule run_raxml:
#     input:
#         phy=twisst_dir + "/gene_vcfs/{sample}.vcf.phylip"
#     output:
#         #best_model= twisst_dir + "/gene_trees/{sample}.raxml.bestModel",
#         best_tree=twisst_dir + "/gene_trees/{sample}.raxml.bestTree",
#         #best_tree_collapsed=twisst_dir + "/gene_trees/{sample}.raxml.bestTreeCollapsed",
#         #bootstraps=twisst_dir + "/gene_trees/{sample}.raxml.bootstraps",
#         #log=twisst_dir + "/gene_trees/{sample}.raxml.log",
#         #ml_trees=twisst_dir + "/gene_trees/{sample}.raxml.mlTrees",
#         #rba=twisst_dir + "/gene_trees/{sample}.raxml.rba",
#         #reduced_phy=twisst_dir + "/gene_trees/{sample}.raxml.reduced.phy",
#         #start_tree=twisst_dir + "/gene_trees/{sample}.raxml.startTree",
#         support=twisst_dir + "/gene_trees/{sample}.raxml.support"
#     threads: 10
#     conda: f"{proj_dir}/envs/scan_phylo.yml"
#     shell:
#         """
#         raxml-ng \
#             --all \
#             --msa {input.phy} \
#             --msa-format PHYLIP \
#             --model GTR \
#             --tree pars{{25}},rand{{25}} \
#             --bs-trees 100 \
#             --threads {threads} \
#             --prefix {twisst_dir}/gene_trees/{wildcards.sample} \
#         """

rule run_raxml:
    input:
        phy=twisst_dir + "/gene_phylips/{sample}.phy"
    output:
        #best_model= twisst_dir + "/gene_trees/{sample}.raxml.bestModel",
        best_tree=twisst_dir + "/gene_trees/{sample}.raxml.bestTree",
        #best_tree_collapsed=twisst_dir + "/gene_trees/{sample}.raxml.bestTreeCollapsed",
        #log=twisst_dir + "/gene_trees/{sample}.raxml.log",
        #ml_trees=twisst_dir + "/gene_trees/{sample}.raxml.mlTrees",
        #rba=twisst_dir + "/gene_trees/{sample}.raxml.rba",
        #reduced_phy=twisst_dir + "/gene_trees/{sample}.raxml.reduced.phy",
        #start_tree=twisst_dir + "/gene_trees/{sample}.raxml.startTree",
    threads: 10
    conda: f"{proj_dir}/envs/scan_phylo.yml"
    shell:
        """
        raxml-ng \
            --search \
            --msa {input.phy} \
            --msa-format PHYLIP \
            --model GTR \
            --tree pars{{50}},rand{{50}} \
            --threads {threads} \
            --prefix {twisst_dir}/gene_trees/{wildcards.sample} \
        """
        
rule newick:
    input:
        gene_trees=expand(twisst_dir + "/gene_trees/{sample}.raxml.bestTree", sample=samples)
    output:
        newick = twisst_dir + "/gene_trees.nwk"
    threads: 1
    conda: f"{proj_dir}/envs/scan_phylo.yml"
    shell:
        """
        cat {input.gene_trees} >{output.newick}
        """
        
# rule collapse:
#     input:
#         newick = twisst_dir + "/gene_trees.nwk"
#     output:
#         newick = twisst_dir + "/gene_trees_BS10.nwk"
#     threads: 1
#     conda: f"{proj_dir}/envs/scan_phylo.yml"
#     shell:
#         """
#         nw_ed {input.newick} 'i & b<=10' o > {output.newick}
#         """
