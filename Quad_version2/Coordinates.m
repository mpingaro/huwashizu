% by Marco Pingaro & Paolo Venini

function [coordinates,nnod]=Coordinates(ndx,ndy,dx,dy)

nnod = (ndx+1)*(ndy+1);
l = ndx*dx; h = ndy*dy;
vecx = 0:dx:l; 
vecy = 0:dy:h;

x = zeros(ndy+1,ndx+1);
y = zeros(ndy+1,ndx+1);

for i=1:ndy+1; x(i,:) = vecx; end
for j=1:ndx+1; y(:,j) = vecy; end

coordinates = zeros(nnod,2);
coordinates(:,1) = reshape(x',nnod,1);
coordinates(:,2) = reshape(y',nnod,1);

return
