# Network data sources

We rely on varying sources for network data.  We are not keeping any
network data in this repository, but we do keep pointers to data sets
and (where appropriate) loader scripts.

## [Koblenz repository][koblenz]

The Koblenz repository includes 253 networks, each of which can be
read from a compressed Tab-Separated Values (TSV) file.  This format
is described in section 9 of the [KONECT handbook][konect-hb].
The bulk of the data is in a file `out.NAME` where `NAME` is the
name of the network; after a couple header lines, this basically
looks like edge information in tuple form.  We provide a reader
for the format under the utilities subdirectory (TODO).

## [UF sparse collection][ufsparse]

The University of Florida Sparse Matrix collection has a wide
variety of test matrices, including many adjacency matrices from
sparse networks.  The easiest way to fetch the matrices is
using the [UFGet](https://www.cise.ufl.edu/research/sparse/mat/UFget.html)
interface from MATLAB.  The `Problem.kind` field identifies graph
or multigraph problems (see [here](https://www.cise.ufl.edu/research/sparse/matrices/kind.html)).
The UF sparse collection contains many of the
[SNAP collection matrices][ufsnap]; when these matrices are
available via UFGet, we recommend it.

## [SNAP repository][snap-data]

The Stanford Large Network Dataset Collection by Jure Leskovec includes
several large data sets.  Unfortunately, they are not in a common format.
Many of these are available
from [the UF sparse collection][ufsnap].

[koblenz]: http://konect.uni-koblenz.de/
[ufsparse]: https://www.cise.ufl.edu/research/sparse/matrices/
[snap-data]: https://snap.stanford.edu/data/
[konect-hb]: http://konect.uni-koblenz.de/downloads/konect-handbook.pdf
[ufsnap]: https://www.cise.ufl.edu/research/sparse/matrices/SNAP/index.html
