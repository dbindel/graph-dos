# UF Sparse Matrix Collection

Assuming the [UFget][ufget] package is installed on your machine,
this loader will fetch remote matrices from the collection.

    % Example use
    [A, problem] = load_graph('ufsparse', 'SNAP/as-caida');

There are several relevant [groups][ufsparse-groups];
see in particular
[Barabasi][Barabasi], [Gleich][Gleich], [Kamvar][Kamvar], 
[Newman][Newman], [Pajek][Pajek], [SNAP][SNAP].

[ufget]: https://www.cise.ufl.edu/research/sparse/mat/UFget.html
[ufsparse-groups: https://www.cise.ufl.edu/research/sparse/matrices/groups.html
[Barabasi]: https://www.cise.ufl.edu/research/sparse/matrices/Barabasi/index.html
[Gleich]: https://www.cise.ufl.edu/research/sparse/matrices/Gleich/index.html
[Kamvar]: https://www.cise.ufl.edu/research/sparse/matrices/Kamvar/index.html
[Newman]: https://www.cise.ufl.edu/research/sparse/matrices/Newman/index.html
[Pajek]: https://www.cise.ufl.edu/research/sparse/matrices/Pajek/index.html
[SNAP]: https://www.cise.ufl.edu/research/sparse/matrices/SNAP/index.html
