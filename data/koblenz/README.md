# Koblenz collection data loaders

Use the `fetch_data` shell script to grab one of the networks for download
from the [Koblenz repository][koblenz]; the data directory will be unpacked
as a subdirectory.

For example, to access the WordNet network:

    # From shell (only needed once)
    ./fetch_data wordnet-words

    % From MATLAB
    A = load_graph('koblenz', 'wordnet-words');

[koblenz]: http://konect.uni-koblenz.de/networks/
