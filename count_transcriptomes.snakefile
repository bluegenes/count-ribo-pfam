"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  diamond_blast.snakefile --use-conda # use -n for dry run
"""

# Download single-copy ribosomal genes

import os
import sys
import itertools

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()


#output files:
#${profile}-c97.faa.clstr.txt
# {sample}-ribo-transcriptome_estimates.txt
single_copy_ribosomal_gene_profiles = "single_copy_ribo.txt"

with open(single_copy_ribosomal_gene_profiles, 'r') as f:
    ribo_profiles = [x.strip() for x in f.readlines()]

rule all:
    input: 
        expand("{sample}_ribomatches/{sample}_ribo_transcriptome_estimates.txt", sample = ["sample"])
  #      expand("{sample}_ribomatches/{sample}_{ribo_profile}-c97.faa.clstr.count", ribo_profile=ribo_profiles, sample=["sample"]),
        #expand("{sample}_ribomatches/{sample}_ribo_transcriptome_estimates.txt", sample = ["sample"]),
        #"sample_ribomatches/sample_ribo_transcriptome_estimates.txt"
#        expand("{samplename}_ribomatches/{samplename}_{ribo_profile}-c97.faa.clstr.txt", ribo_profile=ribo_profiles, samplename= ["sample"])

#  wget http://pfam.xfam.org/family/${profile}/alignment/full -O ${profile}-full.sto
rule download_profile:
    input: lambda wildcards: HTTP.remote(f"http://pfam.xfam.org/family/{wildcards.ribo_profile}/alignment/full", static=True, keep_local=True, allow_redirects=True) 
    output: "pfam_profiles/{ribo_profile}-full.sto" 
    log: "logs/{ribo_profile}.download.log"
    shell: "mv {input} {output} 2> {log}"


#  hmmbuild ${profile}-full.hmm ${profile}-full.sto
rule hmmbuild_profile:
    input:
        "pfam_profiles/{ribo_profile}-full.sto"
    output:
        "pfam_profiles/{ribo_profile}-full.hmm"
    conda:
        "count_tx_env.yml"
    shell:
        """
        hmmbuild {output} {input}
        """

#  hmmpress ${profile}-full.hmm
rule hmmpress_profile:
    input:
        "pfam_profiles/{ribo_profile}-full.hmm"
    output:
        "pfam_profiles/{ribo_profile}-full.hmm.h3f",
        "pfam_profiles/{ribo_profile}-full.hmm.h3i",
        "pfam_profiles/{ribo_profile}-full.hmm.h3m",
        "pfam_profiles/{ribo_profile}-full.hmm.h3p"
    conda:
        "count_tx_env.yml"
    shell:
        """
        hmmpress {input}
        """

#  hmmscan -T 100 --tblout ${profile}-full-tbl.txt --domtblout ${profile}-full-domtbl.txt ${profile}-full.hmm ../example_data/sample_CDSs.faa
rule hmmscan_profile:
    input:
        profile="pfam_profiles/{ribo_profile}-full.hmm.h3f",
        data= 'example_data/{sample}_CDSs.faa', 
    output:
        tblout = "{sample}_ribomatches/{sample}_{ribo_profile}-full-tbl.txt",
        domtblout = "{sample}_ribomatches/{sample}_{ribo_profile}-full-domtbl.txt",
    params:
        prof="pfam_profiles/{ribo_profile}-full.hmm"
    conda:
        "count_tx_env.yml"
    shell:
        """
        hmmscan -T 100 --tblout {output.tblout} --domtblout {output.domtblout} {params.prof} {input.data}
        """

#  cat ${profile}-full-tbl.txt | Rscript -e 'writeLines(noquote(read.table("stdin", stringsAsFactors = F)$V3))' > ${profile}-names.txt
#  python extract-hmmscan-matches.py ${profile}-names.txt ../example_data/sample_CDSs.faa > ${profile}.faa
rule get_fasta_matches:
    input:
        tblout = "{sample}_ribomatches/{sample}_{ribo_profile}-full-tbl.txt",
        data = 'example_data/{sample}_CDSs.faa'
    output:
        "{sample}_ribomatches/{sample}_{ribo_profile}.faa"
    conda:
        "count_tx_env.yml"
    shell:
        """
        python extract-hmmscan-matches.py {input.tblout} --fasta {input.data} > {output}
        """

#cd-hit -i ${profile}.faa -o ${profile}-c97.faa -c .97
#  clstr2txt.pl ${profile}-c97.faa.clstr > ${profile}-c97.faa.clstr.txt
#  grep -v id PF00177-c97.faa.clstr.txt > ${profile}-c97.faa.clstr_nohead.txt
rule cluster_cdhit:
    input:
        "{sample}_ribomatches/{sample}_{ribo_profile}.faa"
    output:
        cluster = "{sample}_ribomatches/{sample}_{ribo_profile}-c97.faa.clstr",
        cluster_txt = "{sample}_ribomatches/{sample}_{ribo_profile}-c97.faa.clstr.txt",
    message:
        "now clustering with cd-hit"
    conda:
        "count_tx_env.yml"
    params:
        cluster_basename = "{sample}_ribomatches/{sample}_{ribo_profile}-c97.faa"
    shell:
        """
        cd-hit -i {input} -o {params.cluster_basename} -c .97
        #if test -f {output}; then
        clstr2txt.pl {output.cluster} > {output.cluster_txt}
        #fi
        """

def aggregate_clusters(w):
   clusters = expand("{sample}_ribomatches/{sample}_{ribo_profile}-c97.faa.clstr.txt", ribo_profile=ribo_profiles, sample=w.sample) 
   return clusters


#num_clusters =$(cat {output} | wc -l)
#echo  -e " {wildcards.ribo_profile} \t {{$num_clusters}} " > {output.cluster_count} 
rule count_clusters:
    input: 
        aggregate_clusters
    output:
        "{sample}_ribomatches/{sample}_ribo_transcriptome_estimates.txt"
    #run:
    #    str_input = " ".join({input})
    #    shell("python count_clusters.py {str_input} > {output}")
    conda:
        "count_tx_env.yml"
    script: 
        "count_clusters.py"

#rule bar_plot:
