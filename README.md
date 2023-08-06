---
output:
  html_document: default
  word_document: default
---
``` r
library(IRanges); library(GenomicRanges); library(CustomBioInfoFuctionsHumanGenome)
```

## 1. Introduction:

This Package contains some customized functions for common operations
used in human genome analysis.

## 2. Availability and Installation

The development version of `CustomBioInfoFuctionsHumanGenome` package is
available at
<https://github.com/vinodsinghjnu/CustomBioInfoFuctionsHumanGenome> and
can be installed as

``` r
# install.packages("devtools")
devtools::install_github("vinodsinghjnu/CustomBioInfoFuctionsHumanGenome",build_vignettes = TRUE )
```

## 3. Functions

### 3.1 addGrInformation

#### **Description**

add assembly information to the genomic range.

#### **Usage**

`addGrInformation(gr.f=gr, assmblyName='hg19')`

#### **Arguments**

-   `gr.f`: A genomic range.
-   `assmblyName`: human genome assembly name i.e., hg19 or hg38.
    Default is `hg19`

#### **Details**

This function will add genomic length and assembly name to given genomic
ranges (Human genome only). It will remove the non-standard chromosomes
from genomic ranges and report the bad genomic ranges for the selected
genome assembly.

#### **Value**

returns the input genomic range along with assembly information.

#### **Examples**

``` r
gr <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),  ranges = IRanges(start = c(1,110,105), end = c(100, 120, 150 )))

outGr=addGrInformation(gr.f=gr, assmblyName='hg19')
```

    ## [1] "bad GRs"
    ## GRanges object with 0 ranges and 0 metadata columns:
    ##    seqnames    ranges strand
    ##       <Rle> <IRanges>  <Rle>
    ##   -------
    ##   seqinfo: 24 sequences from 2 genomes (hg19, NA)

``` r
outGr
```

    ## GRanges object with 3 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]     chr1     1-100      +
    ##   [2]     chr1   110-120      -
    ##   [3]     chr1   105-150      +
    ##   -------
    ##   seqinfo: 24 sequences from 2 genomes (hg19, NA)


### 3.2 pctOverlap_Of_FirstGrToSecondGr

#### **Description**

Percent overlap of a GenomicRange with other

#### **Usage**

`pctOverlap_Of_FirstGrToSecondGr(FirstContext=gr1, SecondContext=gr2)`

#### **Arguments**

-   `FirstContext`: GenomicRange object. (query: overlapped)
-   `SecondContext`: GenomicRange object. (overlapped to)

#### **Details**

Percent overlap of a genomic range with other.

#### **Value**

Percent overlap of first genomic range to second genomic range

#### **Examples**

``` r
gr1 <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),  ranges = IRanges(start = c(1,110,105), end = c(100, 120, 150 )))

gr2 <- GRanges(seqnames = "chr1", strand = c("+", "-", "+", "+"),  ranges = IRanges(start = c(1,12,60, 105), end = c(25, 50, 70, 115 )))

pctOverlap_Of_FirstGrToSecondGr(FirstContext=gr1, SecondContext=gr2)
```

    ## [1] 45.85987

### 3.3 emptyChrGranges

#### **Description**

Create empty chromosome GenomicRange object for a given human genome
assembly

#### **Usage**

`emptyChrGranges(assmblyName='hg19'))`

#### **Arguments**

-   `assmblyName`: hg19 or hg38

#### **Details**

Create empty chromosome GenomicRange object for a given human genome
assembly

#### **Value**

Create empty chromosome GenomicRange object for a given human genome
assembly

Note: Output object is labelled with assembly information.

#### **Examples**

``` r
hg_19_Chr.gr=emptyChrGranges('hg19')
```

    ## [1] "bad GRs"
    ## GRanges object with 0 ranges and 0 metadata columns:
    ##    seqnames    ranges strand
    ##       <Rle> <IRanges>  <Rle>
    ##   -------
    ##   seqinfo: 24 sequences from hg19 genome

``` r
hg_19_Chr.gr
```

    ## GRanges object with 24 ranges and 1 metadata column:
    ##        seqnames      ranges strand |                    Seqs
    ##           <Rle>   <IRanges>  <Rle> |          <DNAStringSet>
    ##    [1]     chr1 1-249250621      * | NNNNNNNNNN...NNNNNNNNNN
    ##    [2]     chr2 1-243199373      * | NNNNNNNNNN...NNNNNNNNNN
    ##    [3]     chr3 1-198022430      * | NNNNNNNNNN...NNNNNNNNNN
    ##    [4]     chr4 1-191154276      * | NNNNNNNNNN...NNNNNNNNNN
    ##    [5]     chr5 1-180915260      * | NNNNNNNNNN...NNNNNNNNNN
    ##    ...      ...         ...    ... .                     ...
    ##   [20]    chr20  1-63025520      * | NNNNNNNNNN...NNNNNNNNNN
    ##   [21]    chr21  1-48129895      * | NNNNNNNNNN...NNNNNNNNNN
    ##   [22]    chr22  1-51304566      * | NNNNNNNNNN...NNNNNNNNNN
    ##   [23]     chrX 1-155270560      * | NNNNNNNNNN...NNNNNNNNNN
    ##   [24]     chrY  1-59373566      * | NNNNNNNNNN...NNNNNNNNNN
    ##   -------
    ##   seqinfo: 24 sequences from hg19 genome

### 3.4 getGbins

#### **Description**

Create GenomicRanges object of given the bin size for human genome.

#### **Usage**

`getGbins(assmblyName='hg19', binSize=1000 )`

#### **Arguments**

-   `assmblyName`: hg19 or hg38
-   `binSize`: size of the genomic block

#### **Details**

Create GenomicRanges object of given the bin size for human genome.

#### **Value**

GenomicRanges object of given the bin size

Note: Output object is labelled with assembly information.

#### **Examples**

``` r
hg_19_Bins.gr=getGbins(assmblyName='hg19', binSize=1000 )
```

    ## [1] "bad GRs"
    ## GRanges object with 0 ranges and 0 metadata columns:
    ##    seqnames    ranges strand
    ##       <Rle> <IRanges>  <Rle>
    ##   -------
    ##   seqinfo: 24 sequences from hg19 genome

``` r
hg_19_Bins.gr
```

    ## GRanges object with 3095679 ranges and 1 metadata column:
    ##             seqnames            ranges strand | CpG_counts
    ##                <Rle>         <IRanges>  <Rle> |  <numeric>
    ##         [1]     chr1            1-1000      * |          0
    ##         [2]     chr1         1001-2000      * |          0
    ##         [3]     chr1         2001-3000      * |          0
    ##         [4]     chr1         3001-4000      * |          0
    ##         [5]     chr1         4001-5000      * |          0
    ##         ...      ...               ...    ... .        ...
    ##   [3095675]     chrY 59369001-59370000      * |          0
    ##   [3095676]     chrY 59370001-59371000      * |          0
    ##   [3095677]     chrY 59371001-59372000      * |          0
    ##   [3095678]     chrY 59372001-59373000      * |          0
    ##   [3095679]     chrY 59373001-59373566      * |          0
    ##   -------
    ##   seqinfo: 24 sequences from hg19 genome

### 3.5 DNASeqsForPattern

#### **Description**

Generate all possible DNA sequences of a [Ambiguous nucleotide
sequence](https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide)

#### **Usage**

`DNASeqsForPattern(pat='NYYN')`

#### **Arguments**

-   `pat`: a [ambiguous nucleotide
    sequence](https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide)

#### **Details**

Generate all possible DNA sequences of a [Ambiguous nucleotide
sequence](https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide)

#### **Value**

A vector of all possible DNA sequences for a given Ambiguous nucleotide
sequence

#### **Examples**

``` r
DNA_seqs=DNASeqsForPattern(pat='NYYN')

DNA_seqs
```

    ##  [1] "ACCA" "CTTC" "GCCG" "TTTT" "ATCA" "CCTC" "GTCG" "TCTT" "ACTA" "CTCC"
    ## [11] "GCTG" "TTCT" "ATTA" "CCCC" "GTTG" "TCCT" "ACCC" "CTTG" "GCCT" "TTTA"
    ## [21] "ATCC" "CCTG" "GTCT" "TCTA" "ACTC" "CTCG" "GCTT" "TTCA" "ATTC" "CCCG"
    ## [31] "GTTT" "TCCA" "ACCG" "CTTT" "GCCA" "TTTC" "ATCG" "CCTT" "GTCA" "TCTC"
    ## [41] "ACTG" "CTCT" "GCTA" "TTCC" "ATTG" "CCCT" "GTTA" "TCCC" "ACCT" "CTTA"
    ## [51] "GCCC" "TTTG" "ATCT" "CCTA" "GTCC" "TCTG" "ACTT" "CTCA" "GCTC" "TTCG"
    ## [61] "ATTT" "CCCA" "GTTC" "TCCG"

### 3.6 createDir_delIfExists

#### **Description**

create a directory and delete if it already exists

#### **Usage**

`createDir_delIfExists(dir='testDir')`

#### **Arguments**

-   `dir`: Name of the dir to be created

#### **Details**

create a directory and delete if it already exists

#### **Value**

Create a directory of given name. (It will delete the directory if it is
already existing there)

#### **Examples**

``` r
createDir_delIfExists(dir='testDir')
```

    ## Old Directory has been deletedNew Directory has been created

``` r
dir.exists('testDir')
```

    ## [1] TRUE

### 3.7 context_oligonucsCounts

#### **Description**

oligo-nucleotide counts in within a genomic context.

#### **Usage**

`context_oligonucsCounts(contextGr=hg_38_gr, oligoType='trinucs', ignore.strand=FALSE, assmblyName='hg38')`

#### **Arguments**

-   `contextGr`: GenomicRange object of the genomic context within which
    oligo-nucleotides has to be counted.
-   `oligoType`: dinucs or trinucs or tetranucs
-   `ignore.strand`: genomic context strand information should be
    considered. Default: FALSE
-   `assmblyName`: human genome assembly name (hg19 or hg38). Default:
    hg19

#### **Details**

oligo-nucleotide counts in within a genomic context.

#### **Value**

a vector of oligo-nucleotide counts

#### **Examples**

``` r
data(hg_38_gr)

oligonucs.Counts=context_oligonucsCounts(contextGr=hg_38_gr, oligoType='trinucs', ignore.strand=FALSE, assmblyName='hg38')

oligonucs.Counts
```

    ##     AAA     AAC     AAG     AAT     ACA     ACC     ACG     ACT     AGA     AGC 
    ## 1325123  507680  697059  857418  698085  401071   87610  563602  766871  484544 
    ##     AGG     AGT     ATA     ATC     ATG     ATT     CAA     CAC     CAG     CAT 
    ##  615364  563602  706620  461672  634677  857418  649772  517039  700902  634677 
    ##     CCA     CCC     CCG     CCT     CGA     CGC     CGG     CGT     CTA     CTC 
    ##  630048  447455   94106  615364   76306   81429   94106   87610  451366  581037 
    ##     CTG     CTT     GAA     GAC     GAG     GAT     GCA     GCC     GCG     GCT 
    ##  700902  697059  685228  326166  581037  461672  496828  405020   81429  484544 
    ##     GGA     GGC     GGG     GGT     GTA     GTC     GTG     GTT     TAA     TAC 
    ##  533465  405020  447455  401071  399481  326166  517039  507680  727231  399481 
    ##     TAG     TAT     TCA     TCC     TCG     TCT     TGA     TGC     TGG     TGT 
    ##  451366  706620  677442  533465   76306  766871  677442  496828  630048  698085 
    ##     TTA     TTC     TTG     TTT 
    ##  727231  685228  649772 1325123
    
    
### 3.8 Gbin_ByCGcnts

#### **Description**

Create a GenomicRange of genomic blocks of user specified CpG counts. (human genome)

#### **Usage**

`Gbin_ByCGcnts(CGs_perBin=100, assmblyName='hg19' )`

#### **Arguments**

-   `assmblyName`: hg19 or hg38 or t2t
-   `CGs_perBin`: CpG counts in a genomic block/bin. (even number)
    oligo-nucleotides has to be counted.
-   `addSeq`: if sequence of the bin is required (Default: FALSE)

#### **Details**

Creates a GenomicRange of genomic blocks of user specified CpG counts. (human genome)

#### **Value**

GenomicRanges object of bins with user specified CpG counts

#### **Examples**

``` r
hg_19_CpGBins.gr=Gbin_ByCGcnts(CGs_perBin=100, assmblyName='hg19' )

hg_19_CpGBins.gr
```

    # GRanges object with 564328 ranges and 2 metadata columns:
    #           seqnames            ranges strand |                     seq CpG_counts
    #              <Rle>         <IRanges>  <Rle> |          <DNAStringSet>  <numeric>
    #       [1]     chr1       10469-10761      * | CGCGGTACCC...CGCGCCGGCG        100
    #       [2]     chr1       10766-11094      * | CGCAGAGAGG...CGTGCACGCG        100
    #       [3]     chr1       11105-12758      * | CGTCACGGTG...AGTGGCGTCG        100
    #       [4]     chr1       12773-15190      * | CGGGGCCGGC...CCCAGCACCG        100
    #       [5]     chr1       15207-17585      * | CGGCTGTTTG...ACACCCCTCG        100
    #       ...      ...               ...    ... .                     ...        ...
    #  [564324]     chrY 59355793-59357713      * | CGACCTGGGC...GAGAGCCACG        100
    #  [564325]     chrY 59357736-59360397      * | CGGATCTCTT...TCACAGCCCG        100
    #  [564326]     chrY 59360409-59361720      * | CGATGGCAGC...CCAACCCCCG        100
    #  [564327]     chrY 59361723-59361953      * | CGTAGGCGTG...GCGCGGCGCG        100
    #  [564328]     chrY 59361962-59362400      * | CGCCTGCGCC...GCGGAAAACG        100
    #  -------
    #  seqinfo: 24 sequences from hg19 genome
    
    
### 3.9 largeVariables

#### **Description**

Memory usage of large variable in workspace

#### **Usage**

`largeVariables(n)`

#### **Arguments**

-   `n`: Number of top memory consuming variables required

#### **Details**

Memory of variable in decreasing order

#### **Value**

DataFrame of top memory variables

#### **Examples**

``` r
largeVariables(n=5)

```


### 3.10 makeTracks_of_grangesList

#### **Description**

make tracks list from GenomicRange lists.

#### **Usage**

`makeTracks_of_grangesList(grlist, location, assmblyName)`

#### **Arguments**

-   `grlist`: GenomicRange object of the genomic context within which
    oligo-nucleotides has to be counted.
-   `location`: location on the chromosome as List. ie., list(chr='chr7',from=26700000, to=26750000) 
-   `assmblyName`: human genome assembly name (hg19 or hg38). Default:
    hg19

#### **Details**

Make tracks from genomic ranges objects list

#### **Value**

list of tracks for input list of genomic ranges objects.

#### **Examples**

``` r
# example 1
data(cpgIslands)

mygrlist=list(gr1=cpgIslands,gr2=cpgIslands)
loc=list(chr='chr7',from=26700000, to=26750000)
tracks=makeTracks_of_grangesList(grlist=mygrlist, location=loc, assmblyName='hg19')
plotTracks(tracks, from = loc$from, to = loc$to)
```


``` r
# example 2, only genomic ranges no annotation
gr1 <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),  ranges = IRanges(start = c(1,110,105), end = c(100, 120, 150 )))

gr2 <- GRanges(seqnames = "chr1", strand = c("+", "-", "+", "+"),  ranges = IRanges(start = c(1,12,60, 105), end = c(25, 50, 70, 115 )))

plotTracks(makeTracks_of_grangesList(list(gr1=gr1,gr2=gr2), if_plain=TRUE), shape="box")
```

