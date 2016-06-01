% [A] = load_koblenz(path)
%
% Returns a (weighted) adjacency matrix from a TSV file
% from the Koblenz network archive.
%
% Input:
%   path: Path to the out.NAME TSV file.
%
% Output:
%   A: Weighted adjacency matrix.
%
function [A,data] = load_koblenz(path)

  % Data set types
  %   sym  - undirected
  %   asym - directed
  %   bip  - bipartite

  % Weights:
  %   unweighted       - only one unweighted edge allowed
  %   positive         - multiple unweighted edges
  %   posweighted      - positively weighted
  %   signed           - any nonzero weight (one edge)
  %   multisigned      - multiple edge version of signed
  %   weighted         - ratings, zero values have no special meaning
  %   multiweighted    - multiple ratings edges
  %   dynamic          - edges can appear or disappear, edges are not weighted
  %   multiposweighted - multiple edge version of posweighted

  % Open file
  [fp,errmsg] = fopen(path, 'rt');
  if fp < 0
    % If it's a directory, try out.whatever
    [path, fname] = fileparts(path);
    path = fullfile(path, fname, ['out.' fname]);
    fp = fopen(path, 'rt');
  end

  % Read first (mandatory) header line
  no_weights = {'unweighted', 'positive'};
  header = fgetl(fp);
  [~, header] = strtok(header);
  [kb_type, header] = strtok(header);
  [kb_wts, header] = strtok(header);
  fprintf('Matrix type: %s\n', kb_type);
  fprintf('Weight type: %s\n', kb_wts);

  % Read second (optional) header line
  header = fgetl(fp);
  [marker, remain] = strtok(header);
  if strcmp(marker, '%')
    counts = sscanf(remain, '%d');
    fprintf('Counts: %d %d %d\n', counts);
    tline = fgetl(fp);
  else
    counts = [0];
    tline = fgetl(fp);
  end

  % Process data lines
  k = 0;
  if any(strcmp(kb_wts, no_weights))
    data = zeros(counts(1), 2);
    while ischar(tline)
      k = k+1;
      data(k,:) = sscanf(tline, '%d%d');
      tline = fgetl(fp);
    end
    data = [data, ones(k,1)];
  else
    data = zeros(counts(1), 3);
    while ischar(tline)
      k = k+1;
      data(k,:) = sscanf(tline, '%d%d%f')
      tline = fgetl(fp);
    end
  end

  fclose(fp);

  % Form the weighted graph
  if length(counts) < 3
    A = sparse(data(:,1), data(:,2), data(:,3));
  else
    A = sparse(data(:,1), data(:,2), data(:,3), counts(2), counts(3));
  end

  % Symmetrize if needed
  if strcmp(kb_type, 'sym')
    d = spdiags(A,0)
    A = spdiags(d,0,A+A');
  end

end
