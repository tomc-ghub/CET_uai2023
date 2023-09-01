function sep = csep(X,Y,Z,cG)
% CSEP Is vertex x d-separated from vertex  given Z in directed (cyclic or acyclic) graph G?
% sep = csep(X, Y, Z, G)  
%
% input:
% - X,Y   : (set of) nodes to separate (typically single vertices)
% - Z     : conditioning set
% - cG    : (possibly) cyclic (bi)directed graph
% output:
% - sep   : 1 = yes: x is d-separated from y given Z in cG

% Note on graph encoding:
% - G = NxN ancestral graph encoded as G(i,j) =
%       0      : not adjacent   i     j
%       1      : tail           i --* j     (with * = {-,>,o}
%       2      : arrowhead      i <-* j
%       3      : -                          (reserved for 'unknown' circle mark 'o')
%       4      : 2-cycle  i <-- j + i --> j (2 arcs, not bidirected edge)
% 
% Approach: 
% 1: expand 2-cycles in G to [i -> u_{ij} -> j] + i <- u_{ji} <- j in H
% 2: do standard d/m-separation on H

DEBUG = ~true;
  %
  if DEBUG,
    strZ = '';
    if ~isempty(Z),
      for z = Z, 
        if isempty(strZ)
          strZ = sprintf('[%i',z); 
        else
          strZ = sprintf('%s,%i',strZ,z); 
        end;
      end;
      strZ = sprintf('%s]',strZ); 
    fprintf('Test csep: %d ][ %d | %s \n',X,Y,strZ); 
    end;
  end;
  
  % no need to search for separating sets of adjacent nodes in G
  if find(cG(X,Y) ~= 0,1), sep = false; return; end;

  % Initialize: 
  M = cG;
  N = length(M);
  M(M == 3) = 1;

  % assume separated until connecting path is found
  sep = true;

  % 2: run standard d-separation
  % construct the graph to update
  G = M;
  % all outgoing arcs with no arrowhead at nodes in Z are eliminated from G
  for z = Z
    Vtmp = find(M(z,:) == 1);
    G(z,Vtmp) = 0;
    G(Vtmp,z) = 0;
  end;  % for s

  % start with all X as noncollider-endpoint
  Vhead = [];
  Vtail = X;
  % while new found (but not yet y) keep trying
  while ~isempty([Vhead,Vtail])
    % reset collector variables for next nodes
    Vnewh = [];
    Vnewt = [];
    % loop over into-nodes 
    for v = Vhead
      % find all nodes (still) reachable from v in G
      % process depending on whether v is selection
      if ~isempty(find(Z == v,1)), 
        % v is selection node: (only arrowheads at v in G)
        v_to = find(G(v,:) ~= 0);
        % check M for edge marks
        v_head = v_to(M(v_to,v) == 2);
        v_tail = v_to(M(v_to,v) == 1);  % or complement
      else
        % v is not selection node .. only tails from v
        v_to = find(G(v,:) == 1);
        % check M for edge marks
        v_head = v_to(M(v_to,v) == 2);
        v_tail = v_to(M(v_to,v) == 1);  % or complement
      end;
      % set traversed edge-half to zero in G and add to nodes found
      G(v,v_to) = 0;
      Vnewh = [Vnewh,v_head];
      Vnewt = [Vnewt,v_tail];
    end;  % for Vhead

    % loop over out-of-nodes 
    for v = Vtail
      % find all nodes (still) reachable from v in G
      if ~isempty(find(Z == v,1)), 
        % non-collider always blocked by selection node
        continue;
      else
        % v is not selection node .. all edge-halves traversable
        v_to = find(G(v,:) ~= 0);
        % check M for edge marks
        v_head = v_to(M(v_to,v) == 2);
        v_tail = v_to(M(v_to,v) == 1);  % or complement
      end;
      % set traversed edge-half to zero in G and add to nodes found
      G(v,v_to) = 0;
      Vnewh = [Vnewh,v_head];
      Vnewt = [Vnewt,v_tail];
    end;  % for Vtail

    % gather next nodes reached at heads/tails
    Vhead = myunique2(Vnewh);
    Vtail = myunique2(Vnewt);    

    % check if we reached y (or set Y) in this step
    if ~isempty(myintersect2([Vhead,Vtail], Y)),
      sep = false;
      return;
    end;
  end; % while
  % if here then not connected
  sep = true;
  % if DEBUG, disp(sprintf('Sep1 = true; Step = %d.',step)); end;
  return;

end  % function msep


% =========================================================================
% local functions to optimize speed 
% =========================================================================
function B = myunique2(A)
% MYUNIQUE Sorts unique elements in set (array) A of positive integers (much faster than built-in unique)
% B = myunique(A)
  if isempty(A)
    B = [];
  else
    bits = zeros(1, max(A));
    bits(A) = 1;
    B = find(bits);
  end
end % myunique2

function C = myintersect2(A,B)
% MYINTERSECT Intersection of two sets of positive integers (much faster than built-in intersect)
% C = myintersect(A,B)
  if isempty(A) || isempty(B), 
    C = []; 
  else
    bits = zeros(1, max([A,B]));
    bits(A) = 1;
    C = B(logical(bits(B)));  
  end;
end  % myintersect2



