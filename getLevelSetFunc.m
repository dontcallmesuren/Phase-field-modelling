function[lsm]=getLevelSetFunc(crack)

global node element
numnode=size(node,1);
cfront=crack(end,:);
seg   = crack(end,1:2) - crack(end-1,1:2);   % tip segment. Like direction vector of crack plane
alpha = atan2(seg(2),seg(1));  % inclination angle
QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
%% Level Set functions
x0  = crack(end-1,1); y0 = crack(end-1,2);
x1  = crack(end,1); y1 = crack(end,2);
t   = 1/norm(seg)*seg;
for i = 1 : numnode
    x = node(i,1);
    y = node(i,2);
    l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0); % = 0 if point on line
    lsm(i,1) = phi/l;            % normal LS
    lsm(i,2) = ([x y]-cfront(1,1:2))*t';  % tangent LS % Angle of (x,y) a.b=|a||b|cos(teta)
end
%% End LSM