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
#NZ_PJOO01000034.1       NZ_JACXCF010000024.1    0.00997914      0       682/1000
#NZ_JABWEE020000068.1    NZ_JAXUFS010000016.1    0.00208065      0       918/1000
#NZ_PKSS02000108.1       NZ_JAXUFS010000016.1    0.00420918      0       844/1000
#NZ_QXXP01000130.1       NZ_MUMV01000094.1       0.0028554       0       890/1000
#NZ_JACGGJ010000012.1    NZ_MUMV01000112.1       0.00399629      0       851/1000

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

cat refseq.nr.fa | head -n 10
#>NZ_JACXCD010000018.1
#TGCGCCAGAATGCGCTGTTGCTTGAGCGAAAACCCGCTGTACCGCCATCGAAACGGGCACATACCGGCAGAGTGGCTGTG
#AAAGAAAGTAATCAGCGATGGTGCTCTGACGGGTTTGAGTTCCGCTGTGATAACGGAGAAAAACTGCGGGTCACGTTCGC
#GCTGGACTGCTGTGACCGTGAGGCACTGCACTGGGCGGTCACAACGGGTGGCTTCGACAGTGAAACAGTACAGGACGTCA
#TGCTGGGAGCAGTGGAACGCCGCTTTGGCAGCGAGCTTCCGGCGTCTCCAGTGGAGTGGCTGACGGATAATGGTTCATGC

rm -fr job
```
# 3 Grouping by MinHash
```
mkdir ~/project/plasmid/grouping
cd ~/project/plasmid/grouping

#get ref sketch
cat ../nr/refseq.nr.fa | mash sketch -k 21 -s 1000 -i -p 8 - -o refseq.nr.k21s1000.msh

#split & get nr plasmid 
faops size ../nr/refseq.nr.fa | cut -f 1 | split -l 1000 -a 3 -d - job/

# get nr sketch
find job -maxdepth 1 -type f -name "[0-9]??" | sort | 
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        faops some ../nr/refseq.nr.fa {} stdout | mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
    '

# calculate dist
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        mash dist -p 6 {}.msh refseq.nr.k21s1000.msh > {}.tsv
    '

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 1 '
        cat {}.tsv
    ' > dist_full.tsv

cat dist_full.tsv | tsv-filter --ff-str-ne 1:2 --le 3:0.05 > connected.tsv
head -n 5 connected.tsv
cat connected.tsv | wc -l

#group
mkdir group
cat connected.tsv | perl -nla -F"\t" -MGraph::Undirected -MPath::Tiny -e '
    BEGIN {
        our $g = Graph::Undirected->new;
    }
    $g->add_edge($F[0], $F[1]);

    END {
        my @rare;
        my $serial = 1;
        my @ccs = $g->connected_components;
        @ccs = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map { [$_, scalar( @{$_} ) ] } @ccs;
        for my $cc ( @ccs ){
            my $count = scalar( @{$cc});
            if ($count < 50){
                push @rare, @{$cc};
            }else{
                path(qq{group/$serial.lst})->spew(map {qq{$_\n}} @{$cc});
                $serial++;
            }
        }
        path(qq{group/00.lst})->spew(map {qq{$_\n}} @rare);

        path(qq{grouped.lst})->spew(map {qq{$_\n}} $g->vertices);
    } 
'
#get non-grouped
faops some -i ../nr/refseq.nr.fa grouped.lst stdout | faops size stdin | cut -f 1 > group/lonely.lst
wc -l group/*

find group -maxdepth 1 -type f -name "[0-9]*.lst" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"

        faops some ../nr/refseq.nr.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
        mash dist -p 6 {}.msh {}.msh > {}.tsv
    '


find group -maxdepth 1 -type f -name "[0-9]*.lst.tsv" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        cat {} |
            tsv-select -f 1-3 |
            Rscript -e '\''
                library(reader);
                library(tidyr);
                library(ape);
                pair_dist <- read.tsv(file("stdin"), col_names=F);
                tmp <- pair_dist %>%
                    pivot_wider(names_form = X2, values_form = X3, values_fill = list(X3 = 1.0))
                tmp <- as.matrix(tmp)
                mat <- tmp[, -1]
                rownames(mat) <-tmp[, 1]

                dist_mat <- as.dist(mat)
                clusters <- hclust(dist_mat, method = "ward.D2")
                tree <- as.phylo(clusters)
                write.tree(phy=tree, file="{.}.tree.nwk")

                group <- cutree(clusters, h=0.2) #k=3
                groups <- as.data.frame(group)
                groups$ids <- rownames(groups)
                rownames(groups) <- NULL
                groups <- groups[order(groups$group), ]
                write_tsv(groups, "{.}.groups.tsv")
            '\''
    '
#.lst.group.tsv
#group   ids
#1       NZ_JRRF01000025.1
#1       NZ_JAJQQO010000031.1
#2       NZ_LQXZ01000009.1

#subgroup
mkdir subgroup
cp group/lonely.lst subgroup/

find group -name "*.groups.tsv" | sort |
    parallel -j 1 -k '
        cat {} | sed -e "1d" | xargs -I[] echo "{/.}_[]"
    ' | sed 's/.lst.groups_/_/' | perl -na -F"\t" -MPath::Tiny -e '
        path(qq{subgroup/$F[0].lst})->append(qq{$F[1]});
    '
#00_1.lst
#NZ_JRRF01000025.1
#NZ_JAJQQO010000031.1

#ignore small subgroups
find subgroup -name "*.lst" | sort |
    parallel -j 1 -k '
        lines=$(cat {} | wc -l)
        if ((lines < 5)); then
            echo -e "{}\t$lines"
            cat {} >> subgroup/lonely.lst
            rm {}
        fi
    '

#append ccs
cat ../nr/connected_components.tsv | parallel -j 1 --colsep "\t" '
    file =$(rg -F -l "{1}" subgroup)
    echo {} | tr "[:blank:]" "\n" >> ${file}
'

#remove duplicates
find subgroup -name "*.lst" | sort |
    parallel -j 1 '
        cat {} | sort | uniq > tmp.lst
        mv tmp.lst {}
    '

wc -l subgroup/* | sort -nr | head -n 100

wc -l subgroup/* | perl -pe 's/^\s+//' | tsv-filter -d" " --le 1:10 | wc -l
#278

wc -1 subgroup/* | perl -pe 's/^\s+//' |
    tsv-filter -d" " --ge 1:50
    tsv-filter -d" " --regex '2:\d+' | sort -nr > next.tsv

#next.tsv
#257 subgroup/00_310.lst
#152 subgroup/00_186.lst

wc -l next.tsv
#41

#rm -fr job
```
# 4 Plasmid:prepare
```
mkdir ~/project/plasmid/genomes
mkdir ~/project/plasmid/taxon

cd ~/project/plasmid/grouping

echo -e "#Serial\tGroup\tCount\tTarget" > ../taxon/group_target.tsv

cat next.tsv | parallel -j 4 -k --line-buffer '
    echo >&2 "==> {}"

    GROUP_NAME = {/.}
    TARGET_NAME = $(head -n 1 {} | perl -pe "s/\.\d+//g")
    SERIAL={#}
    COUNT = $(cat {} | wc -l)
    echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}" >> ../taxon/group_target.tsv

    faops order ../nr/refseq.fa {} stdout | faops filter -s stdin stdout > ../genomes/${GROUP_NAME}.fa
'

cat next.tsv | parallel -j 4 -k --line-buffer '
    echo >&2 "==> {}"

    GROUP_NAME = {/.}
    faops size ../genomes/${GROUP_NAME}.fa > ../taxon/${GROUP_NAME}.sizes
'
#Optional: RepeatMasker
#egaz repeatmasker -p 16 ../genomes/*.fa -o ../genomes/

find genomes -maxdepth 1 -type f -name "*.fa" | sort | parallel -j 4 '
    split-name {} {.}
'

find genomes -maxdepth 2 -mindepth 2 -type f -name "*.fa" | sort | parallel -j 4 '
    mkdir -p {.}
    mv {} {.}
' 
``` 
* preseq
```
cd ~/project/plasmid

cat taxon/group_target.tsv | parallel -j 4 -k --colseq '\t' --line-buffer --no-run-if-empty '
    echo -e "==> Group:[{2}]\tTarget:[{4}]\n"

    for name in $(cat taxon/{2}.sizes | cut -f 1);do
        egaz preseq genomes/{2}/${name}
    done
'
```
* Check outliers of lengths
```
cd ~/project/plasmid

cat taxon/*.sizes | cut -f 1 | wc -l

cat taxon/*.sizes | cut -f 2 | paste -sd+ | bc

cat taxon/group_target.tsv | sed -e "1d" | parallel --colseq '\t' --no-run-if-empty --line-buffer -j 4 -k '
    echo -e "==> Group:[{2}]\tTarget:[{4}]\n"

    median=$(cat taxon/{2}.sizes | datamash median 2)
    mad=$(cat taxon/{2}/sizes | datamash mad 2)
    lower_limit=$(bc <<< " ( ${median} - 2 * ${mad}) / 2 " )

    lines=$(tsv-filter taxon/{2}.sizes --le "2:${lower_limit}" | wc -l)
    if (($lines > 0)); then
        echo >&2 " $lines lines to be filtered " 
        tsv-join taxon/{2}.sizes -e -f <(tsv-filter taxon/{2}.sizes --le "2:${lower_limit}" > taxon/{2}.filtered.sizes
        mv taxon/{2}.filtered.sizes taxon/{2}.sizes
    fi
'
cat taxon/*.sizes | cut -f 1 | wc -l
cat taxon/*.sizes | cut -f 2 | paste -sd+ | bc
```
* Rsync to hpcc
```
rsync -avP ~/project/plasmid wangq@202.119.37.251:data/plasmid
```
# 5 Plasmid:Run
```
cd ~/project/plasmid

cat taxon/group_target.tsv | sed -e "1d" | grep "^53" | parallel --colseq '\t' --no-run-if-empty --line-buffer -j 1 -k '
    echo -e "==> Group:[{2}]\tTarget:[{4}]\n"

    egaz template genomes/{2}/{4} $(cat taxon/{2}.sizes | cut -f 1 | grep -v -x "{4}" | xargs -I[] echo "genomes/{2}/[]"
        -multi -o groups/{2}/ --order --parallel 24 -v

    bsub -q mpi -n 24 -J "{2}-1_pair" "bash groups/{2}/1_pair.sh"
    bsub -w "ended({2}-1_pair)" -q mpi -n 24 -J "{2}-3_multi" "bash groups/{2}/3_multi.sh"
'
    
# others
* MinHash
[MinHash](https://github.com/zhengxy03/language/blob/main/bash/README.md#22-mashminhash)
* Ward分层聚类方法
在每一步合并过程中尽可能减少簇内方差的增加。<br>
簇内平方和（WCSS）：
WCSS是衡量簇内离散程度的指标。对于一个簇，WCSS定义为簇内所有点与簇质心之间距离的平方和。
* MPI
MPI（Message Passing Interface，消息传递接口）是一种标准化和面向高性能计算的通信协议，用于在并行计算环境中的不同进程之间传递消息。它允许多个计算节点（通常是不同的计算机或多核处理器上的不同核心）进行通信和数据交换，以协调工作并解决大规模计算问题。
