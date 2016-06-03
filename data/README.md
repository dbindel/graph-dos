# Data imports

The `load_graph` function in this directory changes to a subdirectory
of the data directory and attempts to call the loader defined there
(which should also be called `load_graph`).  Any arguments to this
`load_graph` are passed through to the loader in the subdirectory.
