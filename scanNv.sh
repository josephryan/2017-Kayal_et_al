#split
cat Nematostella.fa | fasta_formatter -w 0 | split -l 1000  -d - Nvsplit
#run
parallel --jobs 10 '../interproscan.sh -i {} -d . -goterms -f TSV -T {}_temp' ::: ./Nvsplit*
#../interproscan.sh -i ./Nvsplit01  -goterms -f TSV
cat *tsv > allNvGOscan.tsv
#rm Nvsplit*
rm -rf *_temp
