function G = multi_to_simple_graph(M)
% encode (possibly) cyclic graph M as simple graph G (i.e. with one edge
% between nodes)

% 1: initialize
N = size(M,1);
G = M;

% find all edge pairs with both i --> j + i <-- j
[R,C] = find(M == 4);
m = size(R,1);
if (m == 0), return; end;

% extend M with extra entries (all zeros for now)
G(N+m,N+m) = 0;

% 2: process each entry: change x --> y to x --> Gi --> y
for i = 1:m
  x = R(i); y = C(i);
  % modify G
  G(x,N+i) = 1;     % add x --> Gni
  G(N+i,x) = 2;     % 
  G(N+i,y) = 1;     % add Gni --> y
  G(y,N+i) = 2;     % 
  G(x,y) = 0;       % remove original entry from x to y
end;  % for i

% or can we do it in one go?
return;
end % function multi_to_simple_graph