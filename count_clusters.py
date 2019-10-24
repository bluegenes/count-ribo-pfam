import os
import sys
import re
import argparse
import screed
from snakemake.shell import shell

input_files = [snakemake.input] if isinstance(snakemake.input, str) else snakemake.input
input_files_str = [str(i) for i in input_files]

output_file = str(snakemake.output)

def count_clusters(cluster_list, outfile):
    with open(outfile, 'w') as out:
        for clusterfile in cluster_list:
            if not os.path.exists(clusterfile):
                continue
            count = 0
            name = (os.path.basename(clusterfile)).rsplit("-c97.faa.clstr.txt")[0].rsplit('_', 1)[1] # get riboprofile
            with open(clusterfile, 'r') as f:
                for line in f:
                    if not line.startswith("id"):
                        count+=1
            if (count > 0 ):
                out.write(name + '\t' + str(count) + '\n')

count_clusters(input_files_str, output_file)

#if __name__ == '__main__':
#    p = argparse.ArgumentParser()
#    p.add_argument('cluster_txt', nargs="+")
#    args = p.parse_args()
   #assert args.output, "must specify location for output configfile using '-o'"
#    sys.exit(count_clusters(args.cluster_txt))
