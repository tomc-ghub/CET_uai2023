function [G,S,A,M] = mk_random_cg(N,avgK,pCycEdg)
% MK_RANDOM_CG - Create random cyclic directed graph over N variables with
% avg. node degree avgK, and configurable options for the strongly
% connected components (SCC).
% input:
% - N       : nr.of nodes
% - avgK    : average node degree (0 <= avgK < 2(N-1))
% - pCycEdg : probability of [two-cycle, acyclic (with order), cyclic (contra_order)] edges
% output:
% - G       : directed cyclic graph in standard encoding (1-->2 + 4<=>4)
% - S       : NxnSCC matrix of strongly connected components
% - A       : NxN matrix of ancestral relations in G: A(i,j)=1 => i is ancestor of j
% - M       : maximal ancestral graph (cyclic)

global DEBUG;

  % 1 - Initialize  
  % check/process input parameters
  if (nargin < 1), disp('ERROR: no size input N'); return; end; 
  if (nargin < 2), avgK = min((N-1),3.5); end;         % default avg degree 3.5
  if (nargin < 3), pCycEdg = []; end;       % empty = fully random
  % note: two-cycle counts as 2 edges

  G = zeros(N);
  % get nr. of target edges (two cycle x->y + x<-y counts as 2)
  nEdges = round(avgK*N/2);     % so max avgK = (N-1) for cyclic, (N-1)/2 for acyclic

  % sample edges (bit inefficient for very large graphs, but ok for N=1000)
  if ~isempty(pCycEdg)
    % ensure normalize prob
    Z = sum(pCycEdg);
    if Z == 0, return; end;
    p_two = pCycEdg(1)/Z;
    p_acy = pCycEdg(1)/Z;
    % p_cyc = pCycEdg(1)/Z; not needed
    
    % configured to encourage more diversity
    twoE = round(p_two*nEdges/2)*2;    % nr. of edges in two-cycles
    acyE = round(p_acy*nEdges);        % nr. sampled acyclic (with order)
    cycE = nEdges - (twoE+acyE);            % nr. sampled contra order
    
    % step 1: sample two-cycles
    Q = triu(ones(N),1);
    idxG = find(Q > 0);
    if (twoE > 0)
      twoEdges = randsample(size(idxG,1),round(twoE/2),false);
      % set edges in matrix G
      G(idxG(twoEdges)) = 4;    
      G = G + G';
      % update for next step
      Q(idxG(twoEdges)) = 0;
      idxG = find(Q > 0);
    end;
    
    % step 2: sample acyclic (with order) edges
    if (acyE > 0),
      % get allowable indices (update Q from previous step)
      acyEdges = randsample(size(idxG,1),acyE,false);
      % set edges in matrix G
      G(idxG(acyEdges)) = 1;
      G = G';
      G(idxG(acyEdges)) = 2;
      G = G';
      % update Q for next step
      Q(idxG(acyEdges)) = 0;
      idxG = find(Q > 0);
    end;
    
    % step 3: sample against order edges
    if (cycE > 0),
      % get allowable indices (update Q from previous step)
      cycEdges = randsample(size(idxG,1),cycE,false);
      % set edges in matrix G
      G(idxG(cycEdges)) = 2;
      G = G';
      G(idxG(cycEdges)) = 1;
      G = G';
    end;
    
  else
    % fully random (most likely ends up with one big cycle)
    Q = eye(N);
    idxG = find(Q == 0);
    % sample cycE edges without replacement
    Edges = randsample(size(idxG,1),nEdges,false);
    % set corresponding entries to 1
    G(idxG(Edges)) = 1;      % alternatively directly on G and blank diag.
    
    % convert to standard graph encoding
    G((G == 0) & (G' == 1)) = 2;    % arrowheads on arcs
    G((G == 1) & (G' == 1)) = 4;    % 2-cycle
  end;

  % random permute (reorder) rows/columns to avoid any bias
  order = randperm(N);
  G = G(order,order);   % just relabelling the vertices
  
  % get cyclic components+ancestral matrix
  [S,A] = get_scc_an(G);

  % get maximal ancestral graph (MAG)
  M = cg_to_mag(G,S,A);
  
  % done
  return
  
end  % function mk_random_cg










