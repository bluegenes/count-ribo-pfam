

# Download single-copy ribosomal genes

for profile in PF00177 PF00318 PF03946 PF00687 PF00297 PF00573 PF00237 PF00189 PF01929 PF00281 PF00410 PF03719 PF00333 PF00416 PF00411 PF00572 PF00380 PF00203 PF00366 PF00204 PF00521 PF08443 PF00763 PF02882 PF01941 PF01195 PF00625 PF00827 PF00163

do
  wget http://pfam.xfam.org/family/${profile}/alignment/full -O ${profile}-full.sto
  hmmbuild ${profile}-full.hmm ${profile}-full.sto
  hmmpress ${profile}-full.hmm
  hmmscan -T 100 --tblout ${profile}-full-tbl.txt --domtblout ${profile}-full-domtbl.txt ${profile}-full.hmm ../example_data/sample_CDSs.faa
  
  cat ${profile}-full-tbl.txt | Rscript -e 'writeLines(noquote(read.table("stdin", stringsAsFactors = F)$V3))' > ${profile}-names.txt
  
  python extract-hmmscan-matches.py ${profile}-names.txt ../example_data/sample_CDSs.faa > ${profile}.faa
  echo "now clustering with cd-hit"
  cd-hit -i ${profile}.faa -o ${profile}-c97.faa -c .97 
   # convert cluster output to a tsv text file
  clstr2txt.pl ${profile}-c97.faa.clstr > ${profile}-c97.faa.clstr.txt
  grep -v id PF00177-c97.faa.clstr.txt > ${profile}-c97.faa.clstr_nohead.txt
  num_clusters=$(cat ${profile}-c97.faa.clstr_nohead.txt | wc -l)
  echo  -e " ${profile} \t ${num_clusters} " >> transcriptome_estimates.txt
  # bar plot!  

done






