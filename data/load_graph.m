% [A] = load_graph(collection, params)
%
% Loads the (weighted) adjacency for a graph from a directory
% associated with the given collection.  Each collection
% provides its own local version of the load_graph function
% that takes the parameter list (usually to specify a particular
% matrix).
%
function varargout = load_graph(collection, varargin)

  % NB: This does something very strange in Octave --
  %     I get complaints that collection is not defined
  %     even after printing out the value of collection!

  saved_dir = pwd;
  [pathstr,name,ext] = fileparts(mfilename('fullpath'));
  cpath = fullfile(pathstr, collection);
  fname = fullfile(cpath, 'load_graph.m');
  if ~exist(fname, 'file')
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
