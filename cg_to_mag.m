function M = cg_to_mag(G,S,A)
% convert cyclic graph G with SCCs S and ancestral matrix A into maximal
% ancestral graph representation 
% Note: output is a MAG but not quite a proper CMAG yet, as virtual
% v-structures are not marked here (next step in calling routine)
%
% input:
% - G     : NxN graph (directed edges)
% - S     : NxnSCC matrix of nodes strongly connected components (column)
% - A     : NxN ancestor matrix A(i,j) = 1 => i is ancestor of j in G
% output:
% - M     : NxN MAG (maximal ancestral graph, not CMAG)

  DEBUG = ~true;        % toggle if needed to verify 
  % 1 - Initialize  
  % check/process input parameters
  if (nargin < 1), disp('ERROR: no graph input G'); return; end; 
  if (nargin < 3), [S,A] = get_scc_an(G); end;   % 

  N = size(G,1);
  % start from original graph G (standard encoding, 1->2 + 4<=>4)
  M = G;
  
  % 1: induced edges
  % add edges x -> y for all x -> z <- y in G with no edge x-y and y in
  % SCC(z)
  for z = 1:N
    % find all parents of z in H (copy of G with arrowheads
    W = find(G(z,:) >= 2);  % either arrowhead (2) or two-cycle (4)
    % all nodes in the same scc as z
    scc_z = find(S(z,:) == 1,1);        % component column index
    SCC_z = find(S(:,scc_z) == 1);      % nodes in it
    % all possible parents y of z in scc_z
    Y = myintersect(W,SCC_z);
    if ~isempty(Y)
      % loop over all potential x
      for x = W
        % consider all possible y 
        for y = Y
          % if x and y not adjacent 
          if (x ~= y) && (G(x,y) == 0)
            % add arc x -> y (could be overridden by x <- y later if 
            % x in scc(y), but will become x -- y anyway) 
            M(x,y) = 1;
            M(y,x) = 2;
          end;
        end;  % for y
      end;  % for x
    end; % if scc
  end;  % for z
    
  % 2: set all edge marks between nodes in the same SCC to tails
  Q = zeros(size(G));
  % loop over each SCC in S
  for n = 1:size(S,2)
    % mark nodes in component n
    scc_n = find(S(:,n) > 0);
    Q(scc_n,scc_n) = n;
    % set all arrowhead marks in the same SCC in M to tails
    M( (M > 1) & (Q == n) ) = 1;
  end; % for n
  
  return;
end  % function cg_to_mag