#single_chrom_loter.py <chrom_id>

import sys
import os
import numpy as np
import loter.locanc.local_ancestry as lc
import pandas as pd
from pathlib import Path
import pickle

proj_dir="/master/nplatt/sch_hae_scan"
results_dir="{}/results".format(proj_dir)
os.chdir("{}/loter".format(results_dir))

#get chrom
chrom = sys.argv[1]

#get hap ids
hap_ids=[]
with open("hap_ids.csv", 'r') as f:
	hap_ids=f.read().splitlines()

#get variant ids
var_ids = []
with open("h_files/{}_variant_ids.csv".format(chrom), 'r') as f:
	var_ids=f.read().splitlines()

#set up files for output
chr_bovis_hs = np.load("h_files/hs_bovis_{}.npy".format(chrom), mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
chr_haem_hs  = np.load("h_files/hs_haem_{}.npy".format(chrom),  mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
chr_query_hs = np.load("h_files/hs_query_{}.npy".format(chrom), mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
	
#run loter
res_loter = lc.loter_smooth(l_H         = [chr_haem_hs, chr_bovis_hs], 
							h_adm       = chr_query_hs, 
							num_threads = 8) 
	
#process output and save
loter_df = pd.DataFrame(res_loter)
loter_df.index=hap_ids
loter_df.index.name="sample_id"
loter_df.columns=var_ids
loter_df=loter_df.T

try:
   os.remove("csv_files/loter_{}_T.csv".format(chrom))
except OSError:
   pass

loter_df.to_csv("csv_files/loter_{}_T.csv".format(chrom), sep=",")
