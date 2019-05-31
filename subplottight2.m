%function S=subplottight(M,N,dx,dy)
function S=subplottight2(M,N,dx,dy)

if nargin < 3 ||isempty(dx)
    dx = 0.05;
end
if nargin < 4 ||isempty(dy)
    dy = 0.05;
end

if length(dx)==1
    dx = [dx dx dx];
elseif length(dx)==2
    dx = [dx(1) dx(2) dx(2)];
end
if length(dy)==1
    dy = [dy dy dy];
elseif length(dy)==2
    dy = [dy(1) dy(2) dy(2)];
end

w = (1-dx(1)-dx(3)-(N-1)*dx(2))/N;
h = (1-dy(1)-dy(3)-(M-1)*dy(2))/M;

S = zeros(M,N);
for i=1:M
    
    for j=1:N
        S(i,j) = subplot('position',[dx(1)+(j-1)*(w+dx(2)),dy(1)+(i-1)*(h+dy(2)),w,h]);
    end
end
