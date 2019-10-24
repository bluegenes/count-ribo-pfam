"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  diamond_blast.snakefile --use-conda # use -n for dry run
Function: Download single-copy ribosomal gene profiles, use hmmer to map your input data to these, then count matches
"""


import os
import sys
import itertools

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

single_copy_ribosomal_gene_profiles = "single_copy_ribo.txt"

with open(single_copy_ribosomal_gene_profiles, 'r') as f:
    ribo_profiles = [x.strip() for x in f.readlines()]

rule all:
    input: 
        expand("{sample}_ribomatches/{sample}_ribo_transcriptome_estimates.txt", sample = ["sample"])

rule download_profile:
    input: lambda wildcards: HTTP.remote(f"http://pfam.xfam.org/family/{wildcards.ribo_profile}/alignment/full", static=True, keep_local=True, allow_redirects=True) 
    output: "pfam_profiles/{ribo_profile}-full.sto" 
    log: "logs/{ribo_profile}.download.log"
    shell: "mv {input} {output} 2> {log}"


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


rule count_clusters:
    input: 
        aggregate_clusters
    output:
        "{sample}_ribomatches/{sample}_ribo_transcriptome_estimates.txt"
    conda:
        "count_tx_env.yml"
    script: 
        "count_clusters.py"

