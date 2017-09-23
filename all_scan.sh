
#run interproscan on NV sequences from all datasets

../interproscan.sh -i ./OF_AG/AG_nemato_full_seqs.fa -d ./OF_AG/ -goterms -f TSV
../interproscan.sh -i ./OF_AG/OF_nemato_full_seqs.fa -d ./OF_AG/ -goterms -f TSV
../interproscan.sh -i chang_Nemato2.fas -d . -goterms -f TSV
../interproscan.sh -i zapata_Nemato2.fas -d . -goterms -f TSV


#run interproscan on all sequences from the NV genome protein models
#split

cat Nematostella.fa | fasta_formatter -w 0 | split -l 1000  -d - Nvsplit
#run
parallel --jobs 10 '../interproscan.sh -i {} -d . -goterms -f TSV -T {}_temp' ::: ./Nvsplit*
cat *tsv > allNvGOscan.tsv
rm Nvsplit*
rm -rf *_temp
