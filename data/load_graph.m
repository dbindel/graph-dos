% [A] = load_graph(collection, params)
%
% Loads the (weighted) adjacency for a graph from a directory
% associated with the given collection.  Each collection
% provides its own local version of the load_graph function
% that takes the parameter list (usually to specify a particular
% matrix).
%
function varargout = load_graph(collection, varargin)

  saved_dir = pwd;
  [pathstr,name,ext] = fileparts(mfilename('fullpath'));
  if ~exist(fullfile(pathstr, collection, 'load_graph.m'), 'file')
    error(sprintf('Could not find loader for %s', collection));
  end
  try
    cd(fullfile(pathstr, collection));
    [varargout{1:nargout}] = load_graph(varargin{:});
  catch me
    cd(saved_dir);
    rethrow(me);
  end
  cd(saved_dir);
