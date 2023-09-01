function [C,A] = get_scc_an(G)
% combines partitioning (directed) graph G into strongly connected 
% components (SCCs), with obtaining all ancestral relations as matrix
% input:
% - G     : NxN directed cyclic graph (standard or adjacency encoding)
% output:
% - C     : NxnSCC matrix with C(i,j)==1 : node i present in component j
% - A     : NxN matrix    with A(i,j)==1 : node i is ancestor of j
% based on Tarjan's algorithm (linear for sparse graphs, O(|V|+|E|) )
% https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
% ========================================================================
% 
% mini example x -> (y <=> z), with 2 SCCs
% G = [0,1,0;
%      0,0,1;
%      0,1,0];
% SCC: [1],[2,3]
%
% example graph: 8 nodes with 3 strongly connected components (see wikipedia-SCC)
%     1 2   4   6   8
% G  = [0,1,0,0,0,0,0,0;   % 1
%       0,0,1,0,1,1,0,0;   % 
%       0,0,0,1,0,0,1,0;   % 
%       0,0,1,0,0,0,0,1;   % 
%       1,0,0,0,0,1,0,0;   % 
%       0,0,0,0,0,0,1,0;   % 
%       0,0,0,0,0,1,0,0;   % 
%       0,0,0,1,0,0,1,0];  % 
% SCC: [1,2,5], [3,4,8], [6,7]
%
% larger example : 16 nodes with 6 SCCs (see wikipedia-SCC)
%     1 2   4   6   8  10  12  14  16
% G  = [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;   % 1
%       1,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0;   % 
%       0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0;   % 
%       0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0;   % 
%       0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0;   % 
%       0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0;   % 6
%       0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;   % 7
%       0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1;   % 
%       0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1;   % 
%       0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;   % 10
%       0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0;   % 
%       0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;   % 
%       0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0;   % 
%       0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,1;   % 14
%       0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0;   % 
%       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];  % 16
% SCCs: [1,2,3,4,5], [6], [7,8,9], [10,11,12,13], [14,15], [16]

% toggle debug info if needed
DEBUG = ~true;

  % initialize
  N  = size(G,1);           % nr. of variables
  C  = zeros(N,N);          % output C(node,component), max.N components preallocated
  nC = 0;                   % nr. of components
  A  = eye(N,N);            % start with ancestral identity matrix 
  
  Index     = zeros(1,N);   % Index(i) == 0 => index of i undefined
  LowLink   = zeros(1,N);
  onStack   = zeros(1,N);
  % stack of nodes ... preallocate and use index to avoid changing 
  S         = zeros(1,2*N);  % preallocate decent size stack 
  idxS      = 0;

  % loop over all nodes
  idx = 1;            
  for x = 1:N
    if (Index(x) == 0),      % x is not yet defined
      % recursive call to process node
      StrongConnect(x);
    end; 
  end;  % loop over all nodes
  
  % clean up unused entries from SCCs
  C(:,nC+1:N) = [];
  return;
  
  % =========================================================
  % inner function for recursion on successors of or do as loop
  function StrongConnect(v)
  % expand reachable nodes from v
    % process 'strongly_connect(v)'
    Index(v)      = idx;
    LowLink(v)    = idx;
    idx           = idx + 1;
    % stack
    idxS          = idxS + 1;
    S(idxS)       = v;
    onStack(v)    = 1;

    % get all children W of v in G
    W = find(G(v,:) == 1 | G(v,:) == 4);
    % process children w in W
    for w = W
      if Index(w) == 0,      % if w is not yet defined
        % successor w has bot been visited: recurse
        StrongConnect(w);
        LowLink(v) = min([LowLink(v), LowLink(w)]);
      elseif (onStack(w) == 1)
        % successor w has already been visited and is on stack and so
        % part of the current SCC
        LowLink(v) = min([LowLink(v),Index(w)]);
      else
        % if w is not on the stack then it is part of an already handled
        % SCC and should be ignored, except for added descendants
        A(v,:) = A(v,:) | A(w,:);
      end;
    end; % for w

    % if v was the root then pop from stack and create another SCC
    if (LowLink(v) == Index(v))
      % start a new strongly connected component+ track descendants
      
      % = remove nodes from stack until you reach v
      idx_v = find(S == v,1);
      Cnodes = S(idx_v:idxS);
      % add / collect all descendants of nodes in Cnodes
      deC = zeros(1,N);        
      for z = Cnodes
        deC = deC | A(z,:);
      end;
      % assign to all nodes in the component
      for z = Cnodes
        A(z,:) = deC;
      end;
      % SCC is removed from stack, so add descendants to parent entry in stack
      if (idx_v > 1)
        A(S(idx_v-1),:) = A(S(idx_v-1),:) | deC;
      end;
      
      % all nodes in 
      nC = nC + 1;
      % add nodes to new component
      C(Cnodes,nC) = 1;
      % clear up stack
      S(idx_v:idxS) = 0;
      onStack(Cnodes) = 0;
      idxS = idxS - length(Cnodes);
      % report new SCC
      if DEBUG,
        strCnodes = sprintf('%i',x); 
        for c = mysetdiff(Cnodes(:)',x), strCnodes = sprintf('%s,%i',strCnodes,c); end;
        fprintf('Component %d, with Cnodes = [%s]\n',nC, strCnodes); 
      end;
    end;  % root

  end;  % subfunction strongconnect(v)

end % function get_scc_an