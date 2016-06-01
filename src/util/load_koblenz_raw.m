% [A] = load_koblenz_raw(path)
%
% Returns a raw representation of the data in a TSV file
% from the Koblenz network archive.
%
% Input:
%   path: Path to the out.NAME TSV file.
%
% Output:
%   kb: structure with the following fields
%       kb.type   = data set type
%       kb.wts    = weight type
%       kb.counts = triple of counts (relations, src/dest) or empty
%       kb.data   = raw tabular data
%
% The meaning of the data set and weight types is given in Section 9
% of the Koblenz KONECT manual:
%
% Data set types
%   sym  - undirected
%   asym - directed
%   bip  - bipartite
%
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
%
function kb = load_koblenz_raw(path)

  % Open file
  [fp,errmsg] = fopen(path, 'rt');
  if fp < 0
    % If it's a directory, try out.whatever
    [path, fname] = fileparts(path);
    path = fullfile(path, fname, ['out.' fname]);
    fp = fopen(path, 'rt');
  end

  % Read first (mandatory) header line
  header = fgetl(fp);
  [~, header] = strtok(header);
  [kb_type, header] = strtok(header);
  [kb_wts, header] = strtok(header);

  % Read second (optional) header line
  header = fgetl(fp);
  [marker, remain] = strtok(header);
  if strcmp(marker, '%')
    counts = sscanf(remain, '%d');
    hlines = 2;
  else
    counts = [];
    hlines = 1;
  end

  % Wrap up header phase
  fclose(fp);

  % Pack up header data
  kb = [];
  kb.type   = kb_type;
  kb.wts    = kb_wts;
  kb.counts = counts;

  % Read the remainder via dlmread (faster than looping)
  kb.data   = dlmread(path, '', hlines, 0);

end
