## Getting Started

```sh
# Compile the test program for lv89 and edlib
git clone https://github.com/lh3/lv89
cd lv89
make
./ed-test -c test/t2-t.fa test/t2-q.fa

# Compile with WFA2
git clone https://github.com/lh3/lv89
cd lv89
git clone https://github.com/smarco/WFA2-lib
(cd WFA2-lib && make)
make WFA2_ROOT=WFA2-lib
```

## Introduction

This repo provides a reasonably fast yet lightweight implementation of the
[Landau-Vishkin][lv89]/[Myers86][myers86] algorithm to compute the edit
distance between two strings, in a formulation inspired by [WFA2][WFA2]. This
repo doesn't aim at the fastest implementation. It is usually slower than
[edlib][edlib] and [WFA2][WFA2] but is simpler. It also has a unique low-memory
mode to generate base alignment, more effective than the current WFA2.

## The algorithm

The basic LV89 algorithm does global alignment. It takes *O*(*nd*) time and
*O*(*d*) space to compute the distance, where *n* is the length of the longer
sequence and *d* is their edit distance. Deriving the base alignment requires
*O*(*d*<sup>2</sup>) space, which could be significant even when this repo uses
only 2 bits per entry.

To alleviate this issue, we employed a divide-and-conquer approach inspired by
[Hirschberg's algorithm][lin-space]. We split the alignment into ceil(*d*/*b*)
blocks where each block, except the last one, has edit distance exactly *b*.
We then perform base alignment inside each block. The new algorithm takes
*O*(*n*(*d*+*b*)) in time and *O*(*d*\*ceil(*d*/*b*)+*b*<sup>2</sup>) in space. When *b*
equals 1000, the space is reduced by more than two orders of magnitude in
comparison to the basic algorithm.

## Brief benchmark

We only use two pairs of sequences for evaluation. The first pair consists of
the two haplotypes of NA19240 around the C4A/C4B gene. They are 100-150kb in
length with an edit distance of 26k. The second pair consists of GRCh38 and
CHM13 around MHC. They are about 5Mb in length with an edit distance of 123kb.
These sequences can be found [via Zenodo][seq-zenodo].

|Method|CIGAR|CMD option|Time MHC (s)|RAM MHC (MB)|Time C4 (s)|RAM C4 (kB)|
|:-----|:---:|:---------|-----------:|-----------:|----------:|---------:|
|lv89  |N    |          |102         |16          |0.36       |1632|
|edlib |N    |-e        |36          |25          |0.24       |1848|
|WFA2  |N    |-w        |38          |47          |1.52       |4424|
|lv89  |Y    |-c        |156         |3659        |0.68       |22440|
|lv89  |Y    |-cs1000   |151         |191         |0.67       |2652|
|edlib |Y    |-ec       |99          |47          |0.65       |3232|
|WFA2-low|Y  |-wcm1     |146         |1207        |4.50       |109132|
|WFA2-med|Y  |-wcm2     |145         |1310        |4.51       |109132|
|WFA2-high|Y |-wcm3     |74          |58937       |2.51       |3123828|

[edlib][edlib] is the overall winner. It is fast and uses less memory.
[WFA2][WFA2] is slower on the C4 sequences probably because lv89 prunes
diagonals that are unlikely to yield better alignments (this is not a
heuristic). WFA2 also uses more memory to generate base alignment. Nonetheless,
WFA2 is the only library here that supports affine gap penalties, which will be
crucial to many biological applications.

[myers86]: https://link.springer.com/article/10.1007/BF01840446
[lv89]: https://doi.org/10.1016/0196-6774(89)90010-2
[edlib]: https://github.com/Martinsos/edlib
[WFA2]: https://github.com/smarco/WFA2-lib
[lin-space]: https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
[seq-zenodo]: https://zenodo.org/record/6056061
