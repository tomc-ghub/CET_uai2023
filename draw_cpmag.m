function [x, y, h] = draw_cpmag(G, labels, node_t, x, y, varargin)
% DRAW_CPMAG    Draw layout for a MAG / CPAG
% 
% Based to a large extent on 'draw_graph' from FullBNT-1.0.4/GraphViz
% Added capability to draw bi-directed edges and circle marks.
% Note: MAG/CPAG encoding differs from adjacency matrix
  % input:
  % MAG = (maximal) ancestral graph encoded in the form  
  %    (Gij, Gji) = (0,0)  : not adjacent    i     j
  %               = (1,1)  : undirected edge i --- j
  %               = (1,2)  : arrow           i --> j
  %               = (2,1)  : arrow           i <-- j
  %               = (2,2)  : bidirected edge i <-> j
  % CPAG = completed ancestral graph encoded in matrix form with 
  %   G(i,j) = 0  : no mark at i (no edge i-j) i     j
  %          = 1  : tail mark at i (to j)      i --* j    (causal link i=>j?)
  %          = 2  : arrowhead at i (to j)      i <-* j    (i not ancestor of j)
  %          = 3  : circle mark at i (to j)    i o-* j    (non-invariant)  
%
% DRAW_LAYOUT		Draws a layout for a graph 
%
%  [X, Y, H] = DRAW_LAYOUT(ADJ, <LABELS, ISBOX, X, Y>)
%
% Inputs :
%	G : Completed PAG (see above)
%       LABELS : Cell array containing labels <Default : '1':'N'>
%       ISBOX : 1 if node is a box, 0 if oval <Default : zeros>
%       X, Y, : Coordinates of nodes on the unit square <Default : calls make_layout>
%
% Outputs :
%	X, Y : Coordinates of nodes on the unit square
%       H    : Object handles 
%
% Usage Example : [x, y] = draw_layout([0 1;0 0], {'Hidden','Visible'}, [1 0]');
%
% h(i,1) is the text handle - color
% h(i,2) is the circle handle - facecolor
%
% See also MAKE_LAYOUT

% Change History :
% Date		Time		Prog	Note
% 13-Apr-2000	 9:06 PM	ATC	Created under MATLAB 5.3.1.29215a (R11.1)


%G(find(G == 3)) = 1;    % for now: turn circle marks into tails .. (incorrect)

N = size(G,1);
if nargin<2,
  labels = cellstr(int2str((1:N)'));
end

if nargin<3,
  node_t = zeros(N,1);
else
  node_t = node_t(:);
end;
  
axis([0 1 0 1]);
set(gca,'XTick',[],'YTick',[],'box','on','XScale','linear','YScale','linear');
% axis('square');
%colormap(flipud(gray));

if nargin<4
  adj = G;
  adj(find(G(:,:) == 3)) = 1;   % turn circle marks into undirected edges
  adj(find(G(:,:) == 2)) = 1;   % remove arrowheads
  adj(find((G(:,:)+G(:,:)') == 4)) = 1;   % but keep birected edges as undirected for layout 
  [x y] = make_layout(adj);
end;

idx1 = find(node_t==0); h1 = []; wd1=[];
if ~isempty(idx1)
  [h1 wd1] = textoval(x(idx1), y(idx1), labels(idx1), varargin{:});
end;

idx2 = find(node_t~=0); h2 = []; wd2 = [];
if ~isempty(idx2)
  [h2 wd2] = textbox(x(idx2), y(idx2), labels(idx2), varargin{:});
end;

wd = zeros(size(wd1,1)+size(wd2,1),2);
if ~isempty(idx1), wd(idx1, :) = wd1;  end;
if ~isempty(idx2), wd(idx2, :) = wd2; end;

% bug: this code assumes [x y] is the center of each box and oval, which 
% isn't exactly true.
h_edge = [];
for i=1:(N-1),
  for k = (i+1):N
    if (G(i,k) == 0), continue; end;    % if no edge, go to next
    % calculate direction and angle in [x,y] layout
    if x(k)-x(i)==0,
      sign = 1;
      if y(i)>y(k), alpha = -pi/2; else alpha = pi/2; end;
    else
      alpha = atan((y(k)-y(i))/(x(k)-x(i)));
      if x(i)<x(k), sign = 1; else sign = -1; end;
    end;
    % 
    dy1 = sign.*wd(i,2).*sin(alpha);   dx1 = sign.*wd(i,1).*cos(alpha);
    dy2 = sign.*wd(k,2).*sin(alpha);   dx2 = sign.*wd(k,1).*cos(alpha);    
    if (G(i,k)==2) && ((G(k,i)==1) || (G(k,i)==3))           % directed edge i --> k
      h = arrow([x(k)-dx2 y(k)-dy2],[x(i)+dx1 y(i)+dy1],'BaseAngle',30);
      % add small circle o to endpoint
      if (G(k,i) == 3), ellipse((x(k)-dx2), (y(k)-dy2), 0.006, 0.009); end;
    elseif ((G(i,k)==1) || (G(i,k)==3)) && (G(k,i)==2)       % directed edge i <-- k
      h = arrow([x(i)+dx1 y(i)+dy1],[x(k)-dx2 y(k)-dy2],'BaseAngle',30);
      if (G(i,k) == 3), ellipse((x(i)+dx1), (y(i)+dy1), 0.006, 0.009); end;
    elseif (G(i,k)==2) && (G(k,i)==2)       % bidirected edge i <-> k
      xmid = x(k)-dx2 + 0.5*((x(i)+dx1) - (x(k)-dx2));
      ymid = y(k)-dy2 + 0.5*((y(i)+dy1) - (y(k)-dy2));
      h = arrow([xmid ymid],[x(k)-dx2 y(k)-dy2],'BaseAngle',30);
      % note: for some reason a second call draws a black arrow instead of a gray one
      h = arrow([xmid ymid],[x(i)+dx1 y(i)+dy1],'BaseAngle',30);
      
    else	   
      h = line([x(i)+dx1 x(k)-dx2],[y(i)+dy1 y(k)-dy2]);
      % add circle marks 'o' to start/endpoint (if needed)
      if (G(k,i) == 3), ellipse((x(k)-dx2), (y(k)-dy2), 0.006, 0.009); end;
      if (G(i,k) == 3), ellipse((x(i)+dx1), (y(i)+dy1), 0.006, 0.009); end;
    end;
%     if G(k,i)==0, % if directed edge
%       h = arrow([x(i)+dx1 y(i)+dy1],[x(k)-dx2 y(k)-dy2],'BaseAngle',30);
%     else	   
%       h = line([x(i)+dx1 x(k)-dx2],[y(i)+dy1 y(k)-dy2]);
% %      adj(k,i)=-1; % Prevent drawing lines twice ?
%     end;
    h_edge = [h_edge h];
  end;
end;

color.box = 'black';
color.text = color.box;
color.edge = [1 1 1]*3/4;
%color.edge = 'green';
if ~isempty(idx1)
  set(h1(:,1),'Color',color.text)
  set(h1(:,2),'EdgeColor',color.box)
end
if ~isempty(idx2)
  set(h2(:,1),'Color',color.text)
  set(h2(:,2),'EdgeColor',color.box)
end
% set(h_edge,'Color',color.edge)

if nargout>2,
  h = zeros(length(wd),2);
  if ~isempty(idx1),
    h(idx1,:) = h1;
  end;
  if ~isempty(idx2),
    h(idx2,:) = h2;
  end;
end;

%%%%%

function [t, wd] = textoval(x, y, str, varargin)
% TEXTOVAL		Draws an oval around text objects
% 
%  [T, WIDTH] = TEXTOVAL(X, Y, STR)
%  [..] = TEXTOVAL(STR)  % Interactive
% 
% Inputs :
%    X, Y : Coordinates
%    TXT  : Strings
% 
% Outputs :
%    T : Object Handles
%    WIDTH : x and y Width of ovals 
%
% Usage Example : [t] = textoval('Visit to Asia?');
% 
% 
% Note     :
% See also TEXTBOX

% Uses :

% Change History :
% Date		Time		Prog	Note
% 15-Jun-1998	10:36 AM	ATC	Created under MATLAB 5.1.0.421
% 12-Mar-2004   10:00 AM        minka   Changed placement/sizing.
%
% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl 

temp = [];
textProperties = {'BackgroundColor','Color','FontAngle','FontName','FontSize','FontUnits','FontWeight','Rotation'};
varargin = argfilter(varargin,textProperties);

if nargin == 1
  str = x;
end
if ~isa(str,'cell') str=cellstr(str); end;
N = length(str);    
wd = zeros(N,2);
for i=1:N,
  if nargin == 1
    [x, y] = ginput(1);
  end
  tx = text(x(i),y(i),str{i},'HorizontalAlignment','center',varargin{:});
  % minka
  [ptc wx wy] = draw_oval(tx);
  wd(i,:) = [wx wy];
  % draw_oval will paint over the text, so need to redraw it
  delete(tx);
  tx = text(x(i),y(i),str{i},'HorizontalAlignment','center',varargin{:});
  temp = [temp;  tx ptc];
end
if nargout>0, t = temp; end;

%%%%%%%%%


function [ptc, wx, wy] = draw_oval(tx, x, y)
% Draws an oval box around a tex object
sz = get(tx,'Extent');
% minka
wy = 2/3*sz(4);
wx = 2/3*sz(3);
x = sz(1)+sz(3)/2;
y = sz(2)+sz(4)/2;
ptc = ellipse(x, y, wx, wy);
set(ptc, 'FaceColor','w');


%%%%%%%%%%%%%

function [p] = ellipse(x, y, rx, ry, c)
% ELLIPSE		Draws Ellipse shaped patch objects
% 
%  [<P>] = ELLIPSE(X, Y, Rx, Ry, C)
% 
% Inputs :
%    X : N x 1 vector of x coordinates
%    Y : N x 1 vector of y coordinates
%    Rx, Ry : Radii
%    C : Color index
%
% 
% Outputs :
%    P = Handles of Ellipse shaped path objects
% 
% Usage Example : [] = ellipse();
% 
% 
% Note     :
% See also 

% Uses :

% Change History :
% Date		Time		Prog	Note
% 27-May-1998	 9:55 AM	ATC	Created under MATLAB 5.1.0.421

% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl 

if (nargin < 2) error('Usage Example : e = ellipse([0 1],[0 -1],[1 0.5],[2 0.5]); '); end;
if (nargin < 3) rx = 0.1; end;
if (nargin < 4) ry = rx; end;
if (nargin < 5) c = 1; end;

if length(c)==1, c = ones(size(x)).*c; end;
if length(rx)==1, rx = ones(size(x)).*rx; end;
if length(ry)==1, ry = ones(size(x)).*ry; end;
  
n = length(x);
p = zeros(size(x));
t = 0:pi/30:2*pi;
for i=1:n,
	px = rx(i)*cos(t)+x(i);
	py = ry(i)*sin(t)+y(i);
	p(i) = patch(px,py,c(i));
end;

if nargout>0, pp = p; end;

%%%%%

function [t, wd] = textbox(x,y,str,varargin)
% TEXTBOX	Draws A Box around the text 
% 
%  [T, WIDTH] = TEXTBOX(X, Y, STR)
%  [..] = TEXTBOX(STR)
% 
% Inputs :
%    X, Y : Coordinates
%    TXT  : Strings
% 
% Outputs :
%    T : Object Handles
%    WIDTH : x and y Width of boxes 
%% 
% Usage Example : t = textbox({'Ali','Veli','49','50'});
% 
% 
% Note     :
% See also TEXTOVAL

% Uses :

% Change History :
% Date		Time		Prog	Note
% 09-Jun-1998	11:43 AM	ATC	Created under MATLAB 5.1.0.421
% 12-Mar-2004   10:00 AM        minka   Changed placement/sizing.
%
% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl 

temp = [];
textProperties = {'BackgroundColor','Color','FontAngle','FontName','FontSize','FontUnits','FontWeight','Rotation'};
varargin = argfilter(varargin,textProperties);

if nargin == 1
  str = x;
end
if ~isa(str,'cell') str=cellstr(str); end;    
N = length(str);
wd = zeros(N,2);
for i=1:N,
  if nargin == 1
    [x, y] = ginput(1);
  end
  tx = text(x(i),y(i),str{i},'HorizontalAlignment','center',varargin{:});
  % minka
  [ptc wx wy] = draw_box(tx);
  wd(i,:) = [wx wy];
  % draw_box will paint over the text, so need to redraw it
  delete(tx);
  tx = text(x(i),y(i),str{i},'HorizontalAlignment','center',varargin{:});      
  temp = [temp; tx ptc];
end;

if nargout>0, t = temp; end;


function [ptc, wx, wy] = draw_box(tx)
% Draws a box around a text object
sz = get(tx,'Extent');
% minka
wy = 1/2*sz(4);
wx = 1/2*sz(3);
x = sz(1)+sz(3)/2;
y = sz(2)+sz(4)/2;
ptc = patch([x-wx x+wx x+wx x-wx], [y+wy y+wy y-wy y-wy],'w');
set(ptc, 'FaceColor','w');



function args = argfilter(args,keep)
%ARGFILTER  Remove unwanted arguments.
% ARGFILTER(ARGS,KEEP), where ARGS = {'arg1',value1,'arg2',value2,...},
% returns a new argument list where only the arguments named in KEEP are
% retained.  KEEP is a character array or cell array of strings.

% Written by Tom Minka

if ischar(keep)
  keep = cellstr(keep);
end
i = 1;
while i < length(args)
  if ~ismember(args{i},keep)
    args = args(setdiff(1:length(args),[i i+1]));
  else
    i = i + 2;
  end
end
