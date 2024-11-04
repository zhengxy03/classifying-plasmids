# classifying-plasmids
# 1 NCBI RefSeq
```
mkdir -p ~/project/plasmid
cd ~/project/plasmid
mkdir refseq

rsync -avP ftp.ncbi.nlm.nih.gov::refseq/release/plasmid/ RefSeq/
split -b 100M genomic.gbff genomic.gbff_part_
perl ~/project/plasmid/Script/gb_taxon_locus.pl genomic.gff_part_aa > refseq_id_seq. csv
rm genomic.gbff

gzip -dcf RefSeq/plasmid.1.1.genomic.fna.gz | grep "^>" | head -n 5

faops n50 -S -C RefSeq/*.genomic.fna.gz

gzip -dcf RefSeq/*.genomic.fna.gz > RefSeq/plasmid.fa
