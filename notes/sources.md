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
for the format (`load_koblenz`).

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

## [Networkrepository.com][networkrepo]

Like Koblenz, the network data repository has an extensive (500+) set of
contributed network data sets.  The data sets are not in any single
specific format.

## [SNAP repository][snap-data]

The Stanford Large Network Dataset Collection by Jure Leskovec includes
several large data sets.  Unfortunately, they are not in a common format.
Many of these are available from [the UF sparse collection][ufsnap].

## [Gephi data][gephi-data]

The [Gephi](https://gephi.org/) graph visualization project includes
pointers to several data sources.  They are in a variety of formats,
but the formats are standard rather than one-off (and therefore there's
a pretty good chance that a standard converter already exists).

## [RODGER][gleich-rodger]

David Gleich has a Repository of Difficult Graph Experiments and Results
(RODGER) with normalized Laplacian eigenvalues associated with several
graphs.  This is a good source of "ground truth" results for comparison
with the stochastic DoS estimators.

## Other resources

There are a variety of other sites that provide data collections, usually
including much more than network data.  In addition to a ton of one-off
network data sets and small collections associated with a single paper,
there are some portals associated with government data or with data that
is (primarily) produced for citations:

### Government data portals

- <http://dataportals.org/>
- <https://www.data.gov/>
- <https://data.gov.uk/>

### Reproducible research and indexing for data citations

- <http://www.re3data.org/>
- <https://www.datacite.org/>
- <https://figshare.com/>
- <https://datahub.io/>

[koblenz]: http://konect.uni-koblenz.de/
[ufsparse]: https://www.cise.ufl.edu/research/sparse/matrices/
[snap-data]: https://snap.stanford.edu/data/
[konect-hb]: http://konect.uni-koblenz.de/downloads/konect-handbook.pdf
[ufsnap]: https://www.cise.ufl.edu/research/sparse/matrices/SNAP/index.html
[networkrepo]: http://networkrepository.com/
[uci-data]: http://networkdata.ics.uci.edu/
[gephi-data]: https://github.com/gephi/gephi/wiki/Datasets
[gleich-rodger]: https://www.cs.purdue.edu/homes/dgleich/rodger/
