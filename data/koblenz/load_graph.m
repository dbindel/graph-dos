% [A,kb] = load_graph(path)
%
% Returns a (weighted) adjacency matrix from a TSV file
% from the Koblenz network archive.
%
% Input:
%   path: Path to the out.NAME TSV file.
%
% Output:
%   A: Weighted adjacency matrix.
%   kb: Structure of raw output data from loader
%
function [A,kb] = load_graph(path)

  % Load raw data
  if isstruct(path)
    kb = path;
  else
    kb = load_koblenz_raw(path);
  end

  % Unpack the data fields
  data = kb.data;
  i = data(:,1);
  j = data(:,2);
  if size(data,2) == 3
    aij = data(:,3);
  else
    aij = ones(size(i));
  end

  % Form the weighted graph
  counts = kb.counts;
  if length(counts) < 3
    A = sparse(i, j, aij);
  else
    A = sparse(i, j, aij, counts(2), counts(3));
  end

  % Symmetrize if needed
  if strcmp(kb.type, 'sym')
    d = spdiags(A,0);
    A = spdiags(d,0,A+A');
  end

end
