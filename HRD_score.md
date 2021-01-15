# HRD Score

Task: Determine HRD Score

Solution: scarHRD package (<https://github.com/sztup/scarHRD#running-on-mouse-genomes>)

Description:

`scarHRD` is an R package which determines the levels of homologous recombination deficiency (telomeric allelic imbalance, loss off heterozygosity, number of large-scale transitions) based on NGS (WES, WGS) data.

The first genomic scar based homologous recombination deficiency measures were produced using SNP arrays. Since this technology has been largely replaced by next generation sequencing it has become important to develop algorithms that derive the same type of genomic scar-scores from next generation sequencing (WXS, WGS) data. In order to perform this analysis, here we introduce the scarHRD R package and show that using this method the SNP-array based and next generation sequencing based derivation of HRD scores show good correlation.

## Workflow overview

A typical workflow of determining the genomic scar scores for a tumor sample has the following steps:

1. Call allele specific copy number profile on paired normal-tumor BAM files. This step has to be executed before running scarHRD. We recommend using Sequenza (Favero et al. 2015) <http://www.cbs.dtu.dk/biotools/sequenza/> for copy number segmentation, Other tools (e.g. ASCAT (Van Loo et al. 2010)) may also be used in this step.
   This step is time-consuming and compute-intensive.

2. Determine the scar scores with scarHRD R package.
   This step only takes a few minutes.

## Sequenza

Preprocessing of input files
In order to obtain precise mutational and aberration patterns in a tumor sample, Sequenza requires a matched normal sample from the same patient. Typically, the following files are needed to get started with Sequenza:

- A BAM file (or a derived pileup file) from the tumor specimen.
- A BAM file (or a derived pileup file) from the normal specimen.
- A FASTA reference genomic sequence file

The normal and tumor BAM files are processed together to generate a seqz file, which is the required input for the analysis. It is possible to generate a seqz starting from other processed data, such as pileup, or VCF files. The available options are described in the sequenza-utils manual pages.

The sequenza-utils command provides various tools; here we highlight only the basic usage:

Process a FASTA file to produce a GC Wiggle track file:

```console
sequenza−utils gc_wiggle −w 50 --fasta hg19.fa -o hg19.gc50Base.wig.gz
```

Process BAM and Wiggle files to produce a seqz file:

```console
sequenza−utils bam2seqz -n normal.bam -t tumor.bam --fasta hg19.fa \
    -gc hg19.gc50Base.wig.gz -o out.seqz.gz
```

Post-process by binning the original seqz file:

```console
sequenza−utils seqz_binning --seqz out.seqz.gz -w 50 -o out small.seqz.gz
```

The `small.seqz.gz` would be the input to `scarHRD`.

## scarHRD

The scarHRD input may be a detailed segmentation file from Sequenza, in case there is a reliable estimation of ploidy of the tumor sample is known, it should be sumbitted in the ploidy argument of the scarHRD function, otherwise ploidy between 1 and 5.5 will be tested:

```R
a<-read.table("/examples/test1.small.seqz.gz", header=T)
head(a)
```

```R
##   chromosome position base.ref depth.normal depth.tumor depth.ratio    Af
## 1       chr1    12975        N            7          20       2.841 1.000
## 2       chr1    13020        A            8          28       3.500 0.964
## 3       chr1    13026        N           15          43       2.964 1.000
## 4       chr1    13038        T           11          35       3.182 0.971
## 5       chr1    13041        A           11          37       3.364 0.946
## 6       chr1    13077        N           26          65       2.465 1.000
##   Bf zygosity.normal GC.percent good.reads AB.normal AB.tumor tumor.strand
## 1  0             hom         60         51         N        .            0
## 2  0             hom         60         28         A   G0.036         G1.0
## 3  0             hom         59         51         N        .            0
## 4  0             hom         59         35         T   C0.029         C1.0
## 5  0             hom         59         37         A   G0.054         G0.5
## 6  0             hom         62         51         N        .            0
```

Usage example:

```R
library("scarHRD")
scar_score("F:/Documents/scarHRD/examples/test1.small.seqz.gz",reference = "grch38", seqz=TRUE)
```

```R
## Preprocessing started...

## Processing chr1: 18 variant calls; 6290 heterozygous positions; 549112 homozygous positions.
## Processing chr2: 22 variant calls; 4934 heterozygous positions; 394216 homozygous positions.

##
  |=================================================================| 100%
## Preprocessing finished
## Determining HRD-LOH, LST, TAI

##      HRD Telomeric AI LST HRD-sum
## [1,]   1            2   0       3

```

```R
scar_score("F:/Documents/scarHRD/examples/test2.txt",reference = "grch38", seqz=FALSE)
```

```R
## Determining HRD-LOH, LST, TAI

##      HRD Telomeric AI LST HRD-sum
## [1,]  25           35  33      93
```

Genomic scar scores
Loss of Heterozygosity (HRD-LOH)
The HRD-LOH score was described based on investigation in SNP-array-based copy number profiles of ovarian cancer (Abkevich et al. 2012). In this paper the authors showed that the samples with deficient BRCA1, BRCA2 have higher HRD-LOH scores compared to BRCA-intact samples, thus this measurement may be a reliable tool to estimate the sample's homologous recombination capacity.
The definition of a sample's HRD-LOH score is the number of 15 Mb exceeding LOH regions which do not cover the whole chromosome. In the first paper publishing HRD-LOH-score (Abkevich et al., 2012) the authors examine the correlation between HRD-LOH-score and HR deficiency calculated for different LOH region length cut-offs. In that paper the cut-off of 15 Mb approximately in the middle of the interval was arbitrarily selected for further analysis. The authors argue that the rational for this selection rather than selecting the cut-off with the lowest p-value is that the latter cut-off is more sensitive to statistical noise present in the data.
In our manuscript we also investigated if this 15 Mb cutoff is appropriate for WXS-based HRD-LOH score.We followed the same principles as Abkievits et al, thus while there was small difference between the p-values for the different minimum length cutoff values, we chose to use the same, 15 Mb limit as Abkevich et al. We also performed Spearman rank correlation between the SNP-array-based and WXS-based HRD-LOH scores for the different cutoff minimum LOH length cutoff (manuscript, Supplementary Figure S3C). Here the 14 Mb and 15 Mb cutoff-based WXS-HRD-LOH score had the highest correlation with the SNP-based HRD score. (0.700 and 0.695 respectively). This result reassured our choice of using the 15 Mb cutoff like in the SNP-array-based HRD-LOH score.

![alt text](https://github.com/sztup/scarHRD/blob/master/vignettes/hrd-loh.svg)

Figure 1.A Visual representation of the HRD-LOH score on short theoretical chromosomes. Figure 1.B: Calculating HRD-LOH from a biallelic copy-number profile; LOH regions a, and c, would both increase the score by 1, while neither b, or d, would add to its value (b, does not pass the length requirement, and d covers a whole chromosome)

Large Scale Transitions (LST)
The presence of Large Scale Transitions in connection with homologous recombination deficiency was first studied in basal-like breast cancer (Popova et al. 2012). Based on SNP-array derived copy number profiles BRCA1-inactivated cases had showed higher number of large scale transitions.
A large scale transition is defined as a chromosomal break between adjacent regions of at least 10 Mb, with a distance between them not larger than 3Mb.

Figure 2.A: Visual representation of the LST score on short theoretical chromosomes. Figure 2.B: Calculating LST scores from a biallelic copy-number profile; events that are marked with green "marked" signs would increase the score, while events marked with red crosses would not. The grey areas represent the centromeric regions. (From left to right; Chromosome 1: the first event passes the definition of an LST, the second bounded by a shorter than 10 Mb segment from the right, the third is bounded by a segment from the left, which extends to the centromere, the fourth’s gap is greater than 3 Mb. Chromosome 2: The first event is a valid LST, the second and third are not because they are bounded by centromeric segments, and the fourth is a valid LST)

Number of Telomeric Allelic Imbalances
Allelic imbalance (AI) is the unequal contribution of parental allele sequences with or without changes in the overall copy number of the region. Our group have previously found, that the number telomeric AIs is indicative of defective DNA repair in ovarian cancer and triple-negative breast cancer, and that higher number of telomeric AI is associated with better response to cisplatin treatment (Birkbak et al. 2012).
The number of telomeric allelic imbalances is the number AIs that extend to the telomeric end of a chromosome. Figure 3.A: Visual representation of the ntAI on short theoretical chromosomes. Figure 3.B: Illustration of possible telomeric allelic imbalances in an allele specific copy number profile.

```

```
