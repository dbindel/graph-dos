% write_gexf(fname, A, opt)
%
% Write a GEXF file for a graph (https://gephi.org/gexf/format/index.html)
%
% Args:
%   fname: Name of the output file
%   A: Sparse adjacency matrix
%   opt: Optional arguments structure
%     .title: Description of the graph
%     .labels: Labels for the nodes
%     .attr: Cell array of node attribute names (stored as fields in opt)
%
function opt = write_gexf(fname, A, opt)
  if nargin < 3, opt = []; end
  opt = opt_default(A, opt);
  fp = fopen(fname, 'w', 'native', 'UTF-8');
  write_header(fp, opt.title);
  write_attrs(fp, opt);
  write_nodes(fp, length(A), opt);
  write_edges(fp, A);
  write_footer(fp);
  fclose(fp);
end

% -----------------------------------------------------------
% opt = opt_default(A, opt)
%
% Set default options
%
function [opt] = opt_default(A, opt)
  if ~isfield(opt, 'title'), opt.title = 'Exported graph'; end
  if ~isfield(opt, 'attr'),  opt.attr = {};                end
  if ~isfield(opt, 'labels')
    opt.labels = {};
    for k = 1:length(A)
      opt.labels{k} = sprintf('%d', k);
    end
  end
end

% -----------------------------------------------------------
% Write header, footer, attributes, nodes, and edges

function write_attrs(fp, opt)
  if length(opt.attr) == 0
    return;
  end
  fprintf(fp, '        <attributes class="node">\n');
  for k = 1:length(opt.attr)
    fprintf(fp, '          <attribute id="%d" title="%s" type="float"/>\n', ...
            k, opt.attr{k});
  end
  fprintf(fp, '        </attributes>\n');
end

function write_nodes(fp, nnodes, opt)
  fprintf(fp, '        <nodes>\n');
  if length(opt.attr) == 0
    for nid = 1:nnodes
      fprintf(fp, '          <node id="%d" label="%s"/>\n', ...
              nid-1, opt.labels{nid});
    end
  else
    for nid = 1:nnodes
      fprintf(fp, '          <node id="%d" label="%s">\n', ...
              nid-1, opt.labels{nid});
      fprintf(fp, '            <attvalues>\n');
      for k = 1:length(opt.attr)
        a = getfield(opt, opt.attr{k});
        fprintf(fp, '              <attvalue for="%d" value="%g"/>\n', k, a(nid));
      end
      fprintf(fp, '            </attvalues>\n');
      fprintf(fp, '          </node>\n');
    end
  end
  fprintf(fp, '        </nodes>\n');
end

function write_edges(fp, A)
  [i,j,aij] = find(triu(A));
  nedge = length(aij);
  fprintf(fp, '        <edges>\n');
  fprintf(fp, '            <edge id="%d" source="%d" target="%d" />\n', ...
          [0:nedge-1; i'-1; j'-1]);
  fprintf(fp, '        </edges>\n');
end

function write_header(fp, description)
  fmt = ['<?xml version="1.0" encoding="UTF-8"?>\n', ...
         '<gexf xmlns="http://www.gexf.net/1.2draft" version="1.2">\n', ...
         '    <meta lastmodifieddate="%s">\n', ...
         '        <creator>MATLAB exporter</creator>\n', ...
         '        <description>%s</description>\n', ...
         '    </meta>\n', ...
         '    <graph mode="static" defaultedgetype="undirected">\n'];
  fprintf(fp, fmt, datestr(now, 'yyyy-mm-dd'), description);
end

function write_footer(fp)
  footer = ['    </graph>\n', ...
            '</gexf>\n'];
  fprintf(fp, footer);
end
