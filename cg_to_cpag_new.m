function [P,vTriples,uTriples,M,nCounts] = cg_to_cpag_new(G)
% New Cyclic-PAG-from-Graph algorithm (UAI, 2023), incl. corr. to part 2
% (stage f) detailed in the arXiv version of the paper.
% 
% input:
% - G       : directed cyclic graph in standard encoding (1-->2 + 4<=>4)
% output:
% - P       : cyclic PAG (vertices + edges)
% - vTriples: vTriples(i,:) = [a,b,c] = dashed-underlined unshielded collider a -> b <- c
% - uTriples: uTriples(i,:) = [a,b,c] remaining virtual collider triples a -> b + [c]
% - M       : cyclic MAG structure (full CMAG = M + vTriples)
% - nCounts : array with nr. of tests/actions/time spent per stage (6 in total)
  DEBUG = ~true;        % toggle for detailed logging of actions
  
  % Initialize
  nCounts = zeros(3,10);
  N = size(G,1);
  P = [];
  % track underlined (subset of virtual) triples (output)
  vTriples = zeros(N,3); % vTriples(i,:) = [A,..B..,C] (pre-allocate N entries for speed, doubles when needed) 
  uTriples = zeros(N,3); % idem, u-structures 
  idx_vT   = 0;
  idx_uT   = 0;
  vTripM   = zeros(N,N,N); % allow for quick virtual collider triple check
                           % 0=no, 1=virtual v-struct, 2=u-struct

  % PART 1: get CMAG (split in steps a-c to track time spent per step)
  % a: get strongly connected component + ancestors (used in steps b+d)
  tic;
  [SCC,An] = get_scc_an(G);
  nCounts(3,1) = toc;  % elapsed time step a
  
  % b: get MAG (= CMAG without underdashed virtual v-structures)
  tic;
  M = cg_to_mag(G,SCC,An);
  nCounts(3,2) = toc;
  
  % c: v-structures + virtual v-structures (underdashed triples in P) (
  tic;
  % set circle skeleton for corresponding PAG
  P = 3*(M > 0);
  % process all unshielded triples in P
  for b = 1:N
    if DEBUG, fprintf('c: check -> %i <-\n',b); end;
    % get possible a nodes with an arc into b in M
    A = find(M(b,:) == 2);
    C = A;
    for a = A
      % find all other nodes as pair not considered already
      C = mysetdiff(C,a);
      for c = C
        nCounts(1,3) = nCounts(1,3) + 1;
        if (P(a,c) == 0),  % unshielded? (no need to test/check sepset!)
          % 1: orient v-structure a -> b <- c (might get 'underdashed' below)
          P(a,b) = 1; P(b,a) = 2;
          P(c,b) = 1; P(b,c) = 2;
          if DEBUG, fprintf('T2c.1: orient %i -> %i <- %i\n',a,b,c); end;
          nCounts(2,3) = nCounts(2,3) + 1;

          % 2: now check virtual v-structure by comparing to original G
          if (G(a,b) == 0 || G(c,b) == 0),
            % virtual edge ... check it is not a descendant of a common child of a and c in G 
            ch_ac = find(G(:,a) == 2 & G(:,c) == 2);
            % if b is not a descendant of any common children 
            if isempty(ch_ac) || isempty(find(An(ch_ac,b) > 0,1)),
              % virtual!: mark as underlined triple
              idx_vT = idx_vT + 1;
              if (idx_vT > size(vTriples,1)), vTriples(2*end,3) = 0; end;  % double preallocated size 
              vTriples(idx_vT,:) = [min(a,c), b, max(a,c)];
              vTripM(a,b,c) = 1;    % do both to avoid order checks
              vTripM(c,b,a) = 1;
              nCounts(2,4) = nCounts(2,4) + 1;  % nr. of tries = nCounts(2,3)
              if DEBUG, fprintf('T2c.2: underline %i -> [--%i--] <- %i\n',a,b,c); end;
            end;
          end;
        end;
      end; % for c
    end; % for a
  end;  % for b
  nCounts(3,3) = toc;
  % Note: after step c CMAG = [M,vTriples] 
  
  % PART 2: finish construction of CPAG 
  % Note: at this stage [P,vTriples] = skeleton + (virtual) v-structures of P
  % d: remaining edges in u-structures (identifiable arcs into cycles)
  tic;
  % loop over nontrivial SCCs (size > 1)
  idx_ntS = find(sum(SCC,1) > 1);
  for i = idx_ntS
    % get nodes in component
    SCC_i = find(SCC(:,i)' > 0);
    % get parents of SCC in M (arcs, but not from within, so use M)
    paSCC_i = find( sum(M(SCC_i,:) == 2,1) );
    
    % loop over parents
    for a = paSCC_i
      % remaining nonadjacent parents (but not a itself)
      idxC = find( M(a,paSCC_i) == 0);
      C = mysetdiff(paSCC_i(idxC),a);   
      if isempty(C), continue; end;
      
      % get SCC_i \ Adj(a)      (see Lemma 2)
      Adj_a = find(P(a,:) > 0);
      Q = mysetdiff(SCC_i,Adj_a);   % Q = all SCC nodes not adjacent to a

      % get all nodes a o-o {X}
      idxX = P(SCC_i,a) == 3;
      X = SCC_i(idxX);
      
      % loop over all edges a o-o x
      for x = X(:)'
        nCounts(1,5) = nCounts(1,5) + 1;
        % get C_x = subset of C not adjacent to x
        C_x = C(M(x,C) == 0); % possible target parents of SCC_i
        % collect nodes in subgraph
        Qx = myunion(Q,[x,C_x]);       % include x, but a not needed
        
        % get undirected subgraph over Qx
        MQx = (M(Qx,Qx) > 0);        % all undirected

        % get index of x and C in subset Qx
        idx_x = find(Qx == x,1);
        idx_C = find( sum(Qx == C_x',1) > 0);  %index of C_x in Qx
        
        % now get all parent nodes C_x reachable from x in MQx
        visited = zeros(1,size(Qx,2));
        X = idx_x;
        done = false;
        while ~done
          visited(X) = 1;     % ensure each node visited only once
          X = find(sum(MQx(X,:) == 1,1) > 0 & ~visited);  % get adj+not visited
          done = isempty(X);
        end;
        
        % orient a->x and collect all parents nodes visited as uTriple <a,x,C>
        C_visited = Qx(idx_C(find(visited(idx_C))));   % 'find' also allows for []
        if ~isempty(C_visited)
          % u-structure(s) <a,x,.,C_visited> found!
          P(a,x) = 1; P(x,a) = 2;
          if DEBUG, fprintf('T2d: orient %i -> %i \n',a,x); end;
          % store all C_x visited as uTriple <a,x,c>
          for c = C_visited
            % store u-structure
            idx_uT = idx_uT + 1;
            if (idx_uT > size(uTriples,1)), uTriples(2*end,3) = 0; end;  % double preallocated size 
            uTriples(idx_uT,:) = [a, x, c];  % note: always a -> x + [c]
            % also in overall collider triple matrix
            vTripM(a,x,c) = 1;    % do both to avoid order checks
            vTripM(c,x,a) = 1;
            nCounts(2,5) = nCounts(2,5) + 1;  % nr. u-structures found
            if DEBUG, fprintf('T2d: u-triple %i -> %i + [%i]\n',a,x,c); end;
          end;  % for c
        end;  % if (C_visited)
        
      end;  % for x
    end;  % for a
  end; % for i
  % clear up empty entries
  vTriples(idx_vT+1:end,:) = [];  % remove superfluous entries
  uTriples(idx_uT+1:end,:) = [];  % remove superfluous entries
  % elapsed time step d
  nCounts(3,4) = toc; 
  
  % e: ancestorship between adjacent (vTriple) d o-o w (v/uTriple)
  tic;
  % loop over all candidate D to check (in vTriples)
  D = myunion(vTriples(:,2),[]);    % set of unique elements
  % loop over D
  for d = D(:)'
    % find all vTriples <.,d,.>
    idxV = find(vTriples(:,2) == d);
    % find all W to check 
    W = find(P(d,:) == 3);          %  find d o-o w (= not fully oriented)
    % loop over W
    for w = W(:)'
      found = false;
      % loop over all vTriples <.,d,.>
      for i = idxV(:)'
        a = vTriples(i,1);          % vTriple a --> d <-- c 
        c = vTriples(i,3);    
        % try to orient d-w
        nCounts(1,6) = nCounts(1,6) + 1;
        % check for virtual collider triple a --> w <-- c
        if (vTripM(a,w,c) > 0)
          % yes: found one => copy edge from M
          P(d,w) = M(d,w);
          P(w,d) = M(w,d);
          nCounts(2,6) = nCounts(2,6) + 1;
          if DEBUG, 
            if P(d,w) == 1 && P(w,d) == 2
              fprintf('Te.1: orient %i --> %i\n',d,w); 
            elseif P(d,w) == 2 && P(w,d) == 1
              fprintf('Te.2: orient %i <-- %i\n',d,w); 
            else
              fprintf('Te.3: orient %i --- %i\n',d,w); 
            end;
          end; % if debug
        end; % if matching vTrip found
        
        % continue with next edge d-w 
        found = true;
        break;
      end;  % for i
    end;  % for w
  end;  % for d
  nCounts(3,5) = toc;  % elapsed time step e
  
  % f: remaining ancestorship between cycles 
  % corresponds to indirect (v/uTriple) b --* .. --* w <-- d (vTriple) in M
  % Note: for efficiency purposes, we process per edge d o-o w in P, each
  % time finding all possible b in B, and then trying to see if we can
  % reach w from any of them via the set of descendants of B / ancestors of
  % w that are not adjacent to d (nor d itself). If so, then there is an
  % ancestral path from some b in B to w, so also an uncovered path
  % ancestral path, and since the node just before w is not adjacent to d
  % in M (or P) also an uncovered path b --* .. --* w <-- d in M
  tic;
  % again loop over all candidate D (in vTriples) to check (re-use from step e)
  for d = D(:)'
    % find all vTriples <.,d,.>
    idxV = find(vTriples(:,2) == d);
    % collect set of w to consider (edge d o-o w in P, d --> w in M, and w
    % in nontrivial cycle (so |SCC(w)| >= 2) 
    W = find(P(d,:) == 3 & M(:,d)' == 2 & (sum(SCC(:,sum(SCC,1)>=2),2) > 0)' ); 
    
    % loop over candidate w
    for w = W(:)'
      found = false;
      nCounts(1,7) = nCounts(1,7) + 1;
      % process d o-o w in P

      % 1: collect all B : 
      % - find all candidate B with vTrip -> d => v/uTrip
      % - filter out non-ancestors of w
      % - remove adjacent to d (except w)
      B = [];
      % loop over all vTriples <.,d,.>
      for i = idxV(:)'
        a = vTriples(i,1);          % vTriple a --> d <-- c 
        c = vTriples(i,3);    
        % check w is not a descendant of a common child of a and c in M
        % (or easier: directly in G)
        ch_ac = find(G(:,a) == 2 & G(:,c) == 2);  % common children of a and c
        % if w is not a descendant of any common children 
        if isempty(ch_ac) || isempty(find(An(ch_ac,w) > 0,1)),
          % ok: get any (other) virtual collider triples <a,.,c>
          aBc = find(vTripM(a,:,c) > 0);
          B = myunion(B,aBc);
        end;
      end;  % for i (vTriples )
      % filter out all nodes from B that are ancestor of D 
      an_d = find(An(:,d) > 0);
      B = mysetdiff(B, an_d);   
      
      % 2: find all ancestors of w
      an_w = find(An(:,w) > 0);
      % filter out all nodes from B that are not ancestor of w 
      B = myintersect(B,an_w);
      % and not d itself, of course 
      B(B == d) = []; 
      % B-d should be uncovered path, so remove all adj to d in M
      adj_d = find(M(d,:) > 0);
      B = mysetdiff(B, adj_d);
      % check if anything left
      if isempty(B), continue; end;
      
      % 3: get all descendants of B (ensures ancestral path from b to w)
      de_B = find(sum(An(B,:),1) > 0);
      % all candidate nodes on possible uncovered anestral path B-d
      C = myintersect(an_w,de_B);
      % to check if w is on it, (again) remove all other nodes adjacent to d 
      C = mysetdiff(C, adj_d);
      % but not w itself of course
      C = [C,w];  % or 'mysetunion', but this is fine 

      % 4: take subgraph Q of M over nodes C 
      Q = M(C,C);       % note: does not include d (not needed)
      % check if starting from any B we can reach w
      % get index of w and B in C
      idx_w = find(C == w,1);          %index of w in C
      idx_B = find( sum(C == B',1) > 0);    %index of B in C
      
      % 5: get connected
      % note: in contrast to step d this could include directed edges!
      visited = zeros(1,size(Q,2));
      X = idx_B;
      done = false;
      while ~done
        visited(X) = 1;
        X = find(sum(Q(X,:) == 1,1) > 0 & ~visited); % so Q()==1 => tail, so not against an arrowhead
        done = isempty(X);
      end;
          
      % 6: check if we reached w
      if ~isempty(find(visited(idx_w),1))
        % yes: found one => copy edge from M (= d --> w)
        P(d,w) = M(d,w);
        P(w,d) = M(w,d);
        nCounts(2,7) = nCounts(2,7) + 1;
        if DEBUG, 
          if P(d,w) == 1 && P(w,d) == 2
            fprintf('Tf.2a: orient %i --> %i\n',d,w); 
          elseif P(d,w) == 2 && P(w,d) == 1
            fprintf('Tf.2b: orient %i <-- %i\n',d,w); % not triggered
          else
            fprintf('Tf.2c: orient %i --- %i\n',d,w); % idem
          end;
        end; % if debug
        % continue with next w on edge d - w 
      end;  % if found B-w 
      
    end;  % for w
  end; % for d
  nCounts(3,6) = toc;  % elapsed time step f

  % done: clear up and return
  return;
  
end  % cg_to_cpag_new
