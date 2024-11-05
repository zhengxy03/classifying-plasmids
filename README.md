# classifying-plasmids
# 1 NCBI RefSeq
```
#download plasmid refseq data
mkdir -p ~/project/plasmid
cd ~/project/plasmid
mkdir RefSeq

rsync -avP ftp.ncbi.nlm.nih.gov::refseq/release/plasmid/ RefSeq/

#extract taxon and locus information from .gbff file use perl script
gzip -dcf  RefSeq/*.genomic.gbff.gz > genomic.gbff
split -b 100M genomic.gbff genomic.gbff_part_
perl ~/project/plasmid/Script/gb_taxon_locus.pl genomic.gbff_part_aa > refseq_id_seq.csv
rm genomic.gbff

#get and read plasmid genoimc sequence
gzip -dcf RefSeq/plasmid.1.1.genomic.fna.gz | grep "^>" | head -n 5

faops n50 -S -C RefSeq/*.genomic.fna.gz

gzip -dcf RefSeq/*.genomic.fna.gz > RefSeq/plasmid.fa
```
# 2 MinHash to get non-redundant plasmids
```
mkdir nr
cd nr

faops size ../RefSeq/plasmid.fa > refseq.sizes

tsv.filter refseq.sizes --le 2:2000 | wc -l

#extract plasmid sequence > 2000()
faops some ../RefSeq/plasmid.fa <(tsv.filter refseq.sizes --gt 2:2000) refseq.fa

# get ref sketch
cat refseq.fa | mash sketch -k 21 -s 1000 -i -p 8 - -o refseq.plasmid.k21s1000.msh

#split plamsid sequence, get sketch and calculate distance
mkdir job
faops size refseq.fa | cut -f 1 | split -l 1000 -a 3 -d - job/

find job -maxdepth 1 -type f -name "[0-9]??" | sort | parallel -j 4 '
    echo >&2 "==> {}"
    faops some refseq.fa {} stdout | mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
'

find job -maxdepth 1 -type f -name "[0-9]??" | sort | parallel -j 4 '
    echo >&2 "==> {}"
    mash dist {}.msh refseq.plasmid.k21s1000.msh > {}.tsv
'

find job -maxdepth 1 -type f -name "[0-9]??" | sort | parallel -j 4 '
    cat {}.tsv | tsv-filter --ff-str-ne 1:2 --le 3:0.01
' > redundant.tsv

head -n 5 redundant.tsv

cat redundant.tsv | wc -l

cat redundant.tsv | perl -nla -F"\t" -MGraph::Undirected -e '
    BEGIN {
        our $g = Graph::Undirected->new;
    }
    $g->add_edge($F[0], $F[1]);

    END {
        for my $cc ( $g->connected_components ){
            print join qq{\t}, sort @{$cc};
        }
    }
' > connected_components.tsv

cat connected_components.tsv | perl -nla -F"\t" -e '
    printf qq{%s\n}, $_ for @F
' > components.list

wc -l connected_components.tsv components.list

faops some -i refseq.fa components.list stdout >refseq.nr.fa
faops some refseq.fa <(cut -f 1 connected_components.tsv) stdout >> refseq.nr.fa