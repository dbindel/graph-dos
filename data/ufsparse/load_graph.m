% Fetch a problem from the UF Sparse Collection (assumes UFget is installed)
%
function [A, problem] = load_graph(pname)

  problem = UFget(pname);
  A = problem.A;
