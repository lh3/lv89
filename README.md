## Getting Started

```sh
# Compile the test program for lv89 and edlib
git clone https://github.com/lh3/lv89
cd lv89
make

# Compile with WFA2
git clone https://github.com/lh3/lv89
cd lv89
git clone https://github.com/smarco/WFA2-lib
(cd WFA2-lib && make)
make WFA2_ROOT=WFA2-lib
```

## Introduction

This repo implements the [Landau-Vishkin algorithm][lv89] to compute the edit
distance between two strings. This is a fast method for highly similar strings.
[The actual implementation](lv89-full.c) follows a simplified [WFA2][WFA2] formulation
rather than the original formulation. It also learns a performance trick from
WFA. For a pair of ~5Mb MHC sequences with ~123k edits, the implementation here
can find the result in 118 seconds, slower than [edlib][edlib] (83 sec). The
extension mode of lv89 is faster (76 sec). WFA2 gives a wrong result on this
example as of now (2022-03-23).

[lv89]: https://doi.org/10.1016/0196-6774(89)90010-2
[edlib]: https://github.com/Martinsos/edlib
[WFA2]: https://github.com/smarco/WFA2-lib
