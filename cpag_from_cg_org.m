function [P,vTriples,nCounts] = cpag_from_cg_org(G)
% Cyclic PAG-from-Graph algorithm in (Richardson,1996,p18)
% input:
% - G       : directed cyclic graph in standard encoding (1-->2 + 4<=>4)
% output:
% - P       : cyclic PAG (vertices + edges)
% - vTriples: vTriples(i,:) = [a,b,c] = dashed-underlined unshielded collider a -> b <- c
% - nCounts : array with nr. of tests/actions per stage

  % 6 stages, implemented as given, tracks count of operations per stage
  DEBUG = ~true;
  
  % Initialize
  nCounts = zeros(3,7);
  N = size(G,1);
  % separating sets
  Sepset    = cell(N,N);
  SupSepset = cell(N,N,N); % subsequent d-separation/d-connection in cycle
  % track underlined triples? Yes! (even output!)
  vTriples = zeros(N,3); % (i,:) = [A,..B..,C] (pre-allocate N entries for speed, doubles when needed) 
  idxT     = 0;
  
  % obtain extended simple graph H (only single edges between nodes)
  % equivalent to G for d-separation tests
  H = multi_to_simple_graph(G); % adds intermediate nodes at the end
  % helper matrix for checking for parents/children
  Gc = G; Gc(G == 4) = 1;   % encode two-cycle as ancestors/children of eachother 
  
  % get ancestors (needed in steps a+d)
  [SCC,An] = get_scc_an(G);       % note: SCC not needed, but no extra cost
  
  % a: skeleton
  tic;
  % set fully connected circle PAG (all edges x o-o y; P(x,y) == 3)
  P = 3*(ones(N) - eye(N));
  
  % loop over all ordered pairs of vertices
  for x = 1:N
    if DEBUG, fprintf('a: checking x = %i\n',x); end;
    for y = 1:N
      % check d-separation (if still needed)
      if (P(x,y) == 0), continue; end;  % already found a d-sep set
      % construct set S
      Ch_x  = find(Gc(x,:) == 1);                   % note: Gc
      An_xy = find(An(:,x) > 0 | An(:,y) > 0);
      S_xy  = myintersect(Ch_x,An_xy);
      
      % construct set T
      S_xy_x = myunion(S_xy,x);
      Pa_Sx = find( sum(Gc(:,S_xy_x) == 1,2) > 0);  % note: Gc
      Pa_Sx_Sxy = myunion(Pa_Sx,S_xy);
      Ch_y =  find(Gc(y,:)' == 1);
      Ch_xy = myintersect(Ch_x,Ch_y);
      De_chxy = find(sum(An(Ch_xy,:) > 0,1) );
      De_chxy_xy = myunion([x,y]', De_chxy);
      T_xy = mysetdiff(Pa_Sx_Sxy,De_chxy_xy);
 
      % test csep in H (single edge version of G)
      sep = csep(x,y,T_xy,H); % include A? (not needed)
      nCounts(1,1) = nCounts(1,1)+1;        % d-sep tests
      if (sep)
        % remove edge, record sepset
        P(x,y) = 0; P(y,x) = 0;
        Sepset{x,y} = T_xy;
        Sepset{y,x} = T_xy;
        nCounts(2,1) = nCounts(2,1)+1;        % d-sep tests
      end;
    end;
  end;
  nCounts(3,1) = toc;  % elapsed time step a

  % b: v-structures
  tic;
  for b = 1:N
    if DEBUG, fprintf('b: check -> %i <-\n',b); end;
    % get possible a nodes with a link to b in P
    A = find(P(b,:) > 0);
    C = A;
    for a = A
      % find all other nodes as pair not considered already
      C = mysetdiff(C,a);
      for c = C
        if (P(a,c) == 0),
          % unshielded triple: test for v-structure
          nCounts(1,2) = nCounts(1,2) + 1;
          if isempty(find(Sepset{a,c} == b,1))
            % orient v-structure a -> b <- c
            P(a,b) = 1; P(b,a) = 2;
            P(c,b) = 1; P(b,c) = 2;
            if DEBUG, fprintf('b: orient %i -> %i <- %i\n',a,b,c); end;
            nCounts(2,2) = nCounts(2,2) + 1;
          end;
        end;
      end; % for c
    end; % for a
  end;  % for b
  nCounts(3,2) = toc;  % elapsed time step b

  
  % c: orient from non-adjacent
  tic;
  for x = 1:N
    if DEBUG, fprintf('c: x = %i\n',x); end;
    % get possible y nodes with a link to x in P
    Y = find(P(x,:) == 3);      % if not already oriented
    for y = Y
      % find all other nodes A not connected to x or y
      A = find(P(x,:) == 0 & P(y,:) == 0);
      for a = A
        % condition (c): x not in Sepset(a,y)
        Z = Sepset{a,y};
        if isempty(find(Z == x,1))
          sep = csep(a,x,Z,H);
          nCounts(1,3) = nCounts(1,3) + 1;
          if ~sep
            % orient x <- y ...
            P(x,y) = 2; P(y,x) = 1;
            % ... and continue with next edge x o-o y
            if DEBUG, fprintf('c: orient %i <- %i\n',x,y); end;
            nCounts(2,3) = nCounts(2,3) + 1;
            break;
          end;
        end;
      end; % for a
    end; % for y
  end;  % for x
  nCounts(3,3) = toc;  % elapsed time step c
  
  % d: SupSepset (note: could store v-structures at step b)
  tic;
  for b = 1:N
    if DEBUG, fprintf('d: SupSepsets on %i\n',b); end;
    % get possible a nodes with a arc into b in P
    A = find(P(b,:) == 2);
    C = A;
    for a = A
      % find all other nodes as pair not considered already
      C = mysetdiff(C,a);
      for c = C
        % check not adjacent, and not found already
        if (P(a,c) == 0) && isempty(SupSepset{a,b,c}),
          % ok: try construct set Q
          Ch_a  = find(Gc(a,:)' == 1);
          An_abc = find(An(:,a) > 0 | An(:,b) > 0 | An(:,c) > 0);
          Q_abc = myintersect(Ch_a,An_abc);

          % construct set R
          Q_abc_a = myunion(Q_abc,a);
          Pa_Qa = find( sum( Gc(:,Q_abc_a) == 1 ,2) > 0);          
          Pa_Qa_Qabc = myunion(Pa_Qa,Q_abc);
          Ch_c =  find(Gc(c,:)' == 1);
          Ch_ac = myintersect(Ch_a,Ch_c);
          De_Ch_ac = find(An(Ch_ac,:) > 0);
          De_Ch_ac_ac = myunion([a,c]', De_Ch_ac);
          R_abc = mysetdiff(Pa_Qa_Qabc,De_Ch_ac_ac);
          R_abc_b = myunion(R_abc,b);

          % test csep in H (single edge version of G)
          sep = csep(a,c,R_abc_b,H); % include A? no, not useful yet ..
          nCounts(1,4) = nCounts(1,4)+1;        % d-sep tests
          if (sep)
            % SupSepset !
            SupSepset{a,b,c} = R_abc_b;
            SupSepset{c,b,a} = R_abc_b;
            % store underlined triple (sort afterwards
            idxT = idxT + 1;
            if (idxT > size(vTriples,1)), vTriples(2*end,3) = 0; end;  % double preallocated size 
            vTriples(idxT,:) = [min(a,c), b, max(a,c)];
            nCounts(2,4) = nCounts(2,4)+1;        % d-sep tests
          end;  % if (sep)
        end;  % if not adjacent/found already
        
      end; % for c
    end; % for a
  end;  % for b
  nCounts(3,4) = toc;  % elapsed time step d

  
  % e: quadruple of pair of imperfect nonconductors
  tic;
  % loop over all underlined triples 
  for i = 1:idxT
    if DEBUG, fprintf('e: vTriple [%i]\n',i); end;
    % get underlined triple + supsepset
    a = vTriples(i,1);
    b = vTriples(i,2);
    c = vTriples(i,3);
    Z = SupSepset{a,b,c};
    % check existence of nodes adjacent to b with v-structures to a,c
    D = find(P(:,b) == 3 & P(:,a) == 2 & P(:,c) == 2);
    if isempty(D), continue; end;
    for d = D
      nCounts(1,5) = nCounts(1,5) + 1;
      if isempty(find(Z == d,1))
        % orient b --> d
        P(b,d) = 1; P(d,b) = 2;
        if DEBUG, fprintf('e1: orient %i --> %i\n',b,d); end;
        nCounts(2,5) = nCounts(2,5) + 1;
      else
        % orient b *-- d
        P(d,b) = 1;
        if DEBUG, fprintf('e2: orient %i *-- %i\n',b,d); end;
        % so counts = (1,5 - 2,5)
      end;
    end;  % for d
  end;  % for i
  nCounts(3,5) = toc;  % elapsed time step e
  
  % f: extended imperfect nonconductor check
  tic;
  % loop over all underlined triples 
  for i = 1:idxT
    if DEBUG, fprintf('f: vTriple [%i]\n',i); end;
    % get underlined triple + supsepset
    a = vTriples(i,1);
    b = vTriples(i,2);
    c = vTriples(i,3);
    Z = SupSepset{a,b,c};
    % check existence of nodes adjacent to b but not to both a and c
    D = find(P(:,b) == 3 & (P(:,a) == 0 | P(:,c) == 0));
    if isempty(D), continue; end;
    for d = D(:)'
      nCounts(1,6) = nCounts(1,6) + 1;
      Z_d = myunion(Z,d);
      % test csep in G
      sep = csep(a,c,Z_d,H); % note typo in paper (test a,d instead pf a,c)
      if ~sep
        % orient b --> d
        P(b,d) = 1; P(d,b) = 2;
        if DEBUG, fprintf('f: orient %i -> %i\n',b,d); end;
        nCounts(2,6) = nCounts(2,6) + 1;
      end;
    end;  % for d
  end;  % for i
  nCounts(3,6) = toc;  % elapsed time step f
    
  % done: clear up and return
  vTriples(idxT+1:end,:) = [];  % remove superfluous entries
  return;
  
end  % cg_to_cpag_R




